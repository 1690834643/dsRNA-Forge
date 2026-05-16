"""
设计任务生命周期管理
控制层入口：接收 GUI 请求，编排计算引擎，管理数据持久化
"""

import json
from typing import Dict, List, Optional, Callable
from dataclasses import dataclass, asdict

from dsforge.database.manager import DatabaseManager
from dsforge.core.sequence import TranscriptomeIndex, generate_candidates, normalize_sequence, VALID_RNA_BASES
from dsforge.core.scoring.base import evaluate_all_rules
from dsforge.core.dicer import predict_dicer_products, calculate_pool_score
from dsforge.core.explain import explain_result, render_region_map
from dsforge.core.offtarget import OffTargetScreener
from dsforge.core.primers import design_t7_primers
from dsforge.core.redundancy import cluster_redundant_results
from dsforge.core.sgrna import (
    design_genotyping_primers,
    design_sgrna_cloning_oligos,
    scan_sgrna_candidates,
    score_sgrna_offtargets,
)
from dsforge.core.thermodynamics import ThermodynamicsCalculator
from dsforge.core.validation import build_validation_hits


class DesignCancelled(Exception):
    """Raised when a running design task is cancelled by the UI worker."""


@dataclass
class DesignConfig:
    """设计任务配置"""
    mode: str = "long_dsRNA"
    length_min: int = 200
    length_max: int = 500
    gc_min: float = 30.0
    gc_max: float = 52.0
    exclude_poly_n: int = 4
    enabled_rules: List[str] = None
    off_target_levels: Dict = None
    thermodynamics: Dict = None
    n_cores: int = 1
    batch_size: int = 100
    deduplicate_results: bool = True
    cluster_overlap_threshold: float = 0.8
    preset: str = "balanced"
    max_raw_candidates: int = 50000

    def __post_init__(self):
        if self.enabled_rules is None:
            self.enabled_rules = ["consensus", "reynolds", "ui_tei"]
        if self.off_target_levels is None:
            self.off_target_levels = {
                "level_1_16bp": True,
                "level_2_20bp": True,
                "level_3_27bp": False,
                "seed_7nt": True,
            }
        if self.thermodynamics is None:
            self.thermodynamics = {
                "rnaduplex": True,
                "rnacofold": True,
                "rnaup_top_n": 20,
            }


class DesignTask:
    """
    设计任务控制器
    纯 Python，无 Qt 依赖
    """

    def __init__(self, db_manager: Optional[DatabaseManager] = None):
        self.db = db_manager or DatabaseManager()
        self.thermo = ThermodynamicsCalculator()

    def run(
        self,
        transcriptome: TranscriptomeIndex,
        target_seq_id: str,
        config: DesignConfig,
        progress_callback: Optional[Callable[[str, float], None]] = None,
    ) -> Dict:
        """
        执行设计任务（单线程版本，用于 CLI 验证和 GUI 单线程模式）

        Args:
            transcriptome: 转录组索引
            target_seq_id: 目标序列 ID
            config: 设计配置
            progress_callback: 进度回调函数(step, percent)

        Returns:
            {
                'task_id': int,
                'mode': str,
                'results': List[Dict],
                'summary': Dict,
            }
        """
        target_seq = transcriptome.get_sequence(target_seq_id)
        if target_seq is None:
            raise ValueError(f"Target sequence '{target_seq_id}' not found in transcriptome")

        # 创建数据库任务记录
        task_id = self.db.create_task(
            mode=config.mode,
            target_seq_id=target_seq_id,
            target_seq=target_seq,
            params=asdict(config),
        )
        self.db.update_task_status(task_id, "running")

        def _report(step: str, percent: float):
            if progress_callback:
                progress_callback(step, percent)

        try:
            if config.mode == "siRNA":
                results = self._run_sirna_mode(
                    task_id, target_seq_id, target_seq, config, transcriptome, _report
                )
            elif config.mode == "DsiRNA":
                results = self._run_dsiRNA_mode(
                    task_id, target_seq_id, target_seq, config, transcriptome, _report
                )
            elif config.mode == "long_dsRNA":
                results = self._run_long_dsRNA_mode(
                    task_id, target_seq_id, target_seq, config, transcriptome, _report
                )
            elif config.mode == "sgRNA":
                results = self._run_sgrna_mode(
                    task_id, target_seq_id, target_seq, config, transcriptome, _report
                )
            else:
                raise ValueError(f"Unknown mode: {config.mode}")

            self.db.update_task_status(task_id, "completed")
            raw_candidates = getattr(self, "_last_raw_candidate_count", len(results))

            return {
                "task_id": task_id,
                "mode": config.mode,
                "target_seq": target_seq,
                "config": config,
                "results": results,
                "summary": {
                    "total_candidates": len(results),
                    "passed_candidates": sum(1 for r in results if r.get("passed", False)),
                    "raw_candidates": raw_candidates,
                    "nonredundant_candidates": len(results),
                    "deduplicated_candidates": max(0, raw_candidates - len(results)),
                },
            }

        except DesignCancelled:
            self.db.update_task_status(task_id, "cancelled")
            raise
        except Exception as e:
            try:
                self.db.update_task_status(task_id, "failed")
            except Exception as db_err:
                raise RuntimeError(
                    f"Design failed ({e}) and additionally database status update failed: {db_err}"
                ) from e
            raise

    def _candidate_window_count(self, target_seq, min_len: int, max_len: int) -> int:
        """Return raw clean sliding windows before GC/repeat filtering."""
        if isinstance(target_seq, int):
            seq = "A" * max(0, target_seq)
        else:
            seq = normalize_sequence(target_seq or "")
        if len(seq) <= 0 or min_len <= 0 or max_len < min_len:
            return 0
        total = 0
        for length in range(min_len, max_len + 1):
            if len(seq) < length:
                continue
            for start in range(0, len(seq) - length + 1):
                window = seq[start : start + length]
                if set(window) <= VALID_RNA_BASES:
                    total += 1
        return total

    def _ensure_candidate_budget(self, target_seq: str, min_len: int, max_len: int, config: DesignConfig):
        max_raw = int(getattr(config, "max_raw_candidates", 0) or 0)
        if max_raw <= 0:
            return
        raw_count = self._candidate_window_count(target_seq or "", min_len, max_len)
        if raw_count > max_raw:
            raise ValueError(
                f"Too many candidate windows ({raw_count:,}) for this target and length range. "
                f"Use a shorter target/region, narrow the length range, or raise max_raw_candidates."
            )

    def _run_sirna_mode(
        self,
        task_id: int,
        target_seq_id: str,
        target_seq: str,
        config: DesignConfig,
        transcriptome: TranscriptomeIndex,
        report: Callable,
    ) -> List[Dict]:
        """siRNA 模式：21nt 候选"""
        report("Generating siRNA candidates...", 0)
        self._ensure_candidate_budget(target_seq, 21, 21, config)

        candidates = list(generate_candidates(
            target_seq,
            mode="siRNA",
            min_len=21,
            max_len=21,
            gc_min=config.gc_min,
            gc_max=config.gc_max,
            exclude_poly=config.exclude_poly_n,
            max_candidates=config.max_raw_candidates,
        ))

        report(f"Generated {len(candidates)} candidates", 10)

        report("Building or reusing off-target risk index...", 12)
        screener = OffTargetScreener(transcriptome)
        screener.prepare_index()
        results = []

        for i, cand in enumerate(candidates):
            seq = cand["sequence"]

            # 多规则评分
            rule_results = evaluate_all_rules(seq, config.enabled_rules)
            consensus = rule_results.get("consensus", {})

            # 脱靶筛查
            off_target = screener.screen_sequence(
                seq, exclude_ids={target_seq_id}, **config.off_target_levels
            )

            # 热力学（简化版：仅 on-target RNAduplex）
            thermo = None
            if self.thermo.available and config.thermodynamics.get("rnaduplex"):
                thermo = self.thermo.rnaduplex(seq, target_seq)

            result_record = {
                "task_id": task_id,
                "rank": i + 1,
                "sequence": seq,
                "position": f"{cand['start']}-{cand['end']}",
                "position_start": cand["start"],
                "position_end": cand["end"],
                "consensus_score": consensus.get("score", 0),
                "passed": consensus.get("passed", False) and off_target["passed"],
                "rules": rule_results,
                "off_target": off_target,
                "thermodynamics": thermo,
            }
            results.append(result_record)

            percent = 10 + (i + 1) / len(candidates) * 70
            report(f"Scoring candidate {i+1}/{len(candidates)}", percent)

        report("Ranking candidates by score and off-target risk...", 84)
        results = self._rank_and_deduplicate(results, config)
        self._refine_top_results_with_rnaup(results, target_seq, transcriptome, config, report, 86, 93)
        self._annotate_decision_artifacts(results, target_seq_id, target_seq, transcriptome, config)
        results.sort(key=lambda x: x.get("recommendation_score", x["consensus_score"]), reverse=True)
        for i, r in enumerate(results):
            r["rank"] = i + 1

        # 保存到数据库
        report("Saving results...", 94)
        self._save_results(task_id, results, config)
        report("Design computation complete", 98)

        return results

    def _rank_and_deduplicate(self, results: List[Dict], config: DesignConfig) -> List[Dict]:
        """Sort and optionally collapse highly overlapping sliding-window neighbors."""
        self._annotate_recommendation_scores(results)
        results.sort(key=lambda x: x.get("recommendation_score", x["consensus_score"]), reverse=True)
        for i, r in enumerate(results):
            r["rank"] = i + 1
            r.setdefault("cluster_id", i + 1)
            r.setdefault("cluster_size", 1)
            r.setdefault("alternative_count", 0)
            r.setdefault("cluster_span", r.get("position", ""))
        self._last_raw_candidate_count = len(results)
        if config.deduplicate_results:
            results = cluster_redundant_results(results, config.cluster_overlap_threshold)
        return results

    def _run_dsiRNA_mode(
        self,
        task_id: int,
        target_seq_id: str,
        target_seq: str,
        config: DesignConfig,
        transcriptome: TranscriptomeIndex,
        report: Callable,
    ) -> List[Dict]:
        """DsiRNA 模式：27nt 候选，预测切割产物"""
        report("Generating DsiRNA candidates...", 0)
        self._ensure_candidate_budget(target_seq, 27, 27, config)

        candidates = list(generate_candidates(
            target_seq,
            mode="DsiRNA",
            min_len=27,
            max_len=27,
            gc_min=config.gc_min,
            gc_max=config.gc_max,
            exclude_poly=config.exclude_poly_n,
            max_candidates=config.max_raw_candidates,
        ))

        report(f"Generated {len(candidates)} candidates", 10)

        report("Building or reusing off-target risk index...", 12)
        screener = OffTargetScreener(transcriptome)
        screener.prepare_index()
        results = []

        for i, cand in enumerate(candidates):
            seq = cand["sequence"]

            # 预测 Dicer 切割产物
            products = predict_dicer_products(seq, cut_length=21)

            # 评估切割产物质量
            pool_result = calculate_pool_score(
                products,
                scorer=evaluate_all_rules,
                enabled_rules=config.enabled_rules,
            )

            # 脱靶筛查（对 pool 中的产物进行）；seed 开关只控制 seed 层级，
            # 不能关闭 16/20/27bp 连续匹配层级。
            off_target = screener.screen_pool(
                products, exclude_ids={target_seq_id}, **config.off_target_levels
            )

            result_record = {
                "task_id": task_id,
                "rank": i + 1,
                "sequence": seq,
                "position": f"{cand['start']}-{cand['end']}",
                "position_start": cand["start"],
                "position_end": cand["end"],
                "consensus_score": pool_result["pool_score"],
                "passed": pool_result["pool_score"] > 50 and off_target["passed"],
                "pool": pool_result,
                "off_target": off_target,
            }
            results.append(result_record)

            percent = 10 + (i + 1) / len(candidates) * 70
            report(f"Scoring DsiRNA {i+1}/{len(candidates)}", percent)

        report("Ranking candidates by score and off-target risk...", 84)
        results = self._rank_and_deduplicate(results, config)
        self._refine_top_results_with_rnaup(results, target_seq, transcriptome, config, report, 86, 93)
        self._annotate_decision_artifacts(results, target_seq_id, target_seq, transcriptome, config)
        results.sort(key=lambda x: x.get("recommendation_score", x["consensus_score"]), reverse=True)
        for i, r in enumerate(results):
            r["rank"] = i + 1

        report("Saving results...", 94)
        self._save_results(task_id, results, config)
        report("Design computation complete", 98)
        return results

    def _run_long_dsRNA_mode(
        self,
        task_id: int,
        target_seq_id: str,
        target_seq: str,
        config: DesignConfig,
        transcriptome: TranscriptomeIndex,
        report: Callable,
    ) -> List[Dict]:
        """长 dsRNA 模式：200-500bp 区段"""
        report("Generating long dsRNA candidates...", 0)
        self._ensure_candidate_budget(target_seq, config.length_min, config.length_max, config)

        candidates = list(generate_candidates(
            target_seq,
            mode="long_dsRNA",
            min_len=config.length_min,
            max_len=config.length_max,
            gc_min=config.gc_min,
            gc_max=config.gc_max,
            exclude_poly=config.exclude_poly_n,
            max_candidates=config.max_raw_candidates,
        ))

        report(f"Generated {len(candidates)} candidate regions", 5)

        report("Building or reusing off-target risk index...", 8)
        screener = OffTargetScreener(transcriptome)
        screener.prepare_index()
        results = []

        for i, cand in enumerate(candidates):
            seq = cand["sequence"]

            # 预测 Dicer 切割产物
            products = predict_dicer_products(seq)

            # Pool 整体评分
            pool_result = calculate_pool_score(
                products,
                scorer=evaluate_all_rules,
                enabled_rules=config.enabled_rules,
            )

            # Pool 脱靶筛查
            pool_offtarget = screener.screen_pool(
                products, exclude_ids={target_seq_id}, **config.off_target_levels
            )

            # 热力学（简化：仅对 pool 中 top 产物计算 on-target）
            thermo = None
            if self.thermo.available and products:
                top_product = max(products, key=lambda p: len(p["sequence"]))
                thermo = self.thermo.rnaduplex(top_product["sequence"], target_seq)

            result_record = {
                "task_id": task_id,
                "rank": i + 1,
                "sequence": seq,
                "position": f"{cand['start']}-{cand['end']}",
                "position_start": cand["start"],
                "position_end": cand["end"],
                "consensus_score": pool_result["pool_score"],
                "passed": pool_result["pool_score"] > 50 and pool_offtarget["pool_risk"] != "high",
                "pool": pool_result,
                "off_target": pool_offtarget,
                "thermodynamics": thermo,
            }
            results.append(result_record)

            percent = 5 + (i + 1) / len(candidates) * 75
            report(f"Evaluating region {i+1}/{len(candidates)}", percent)

        report("Ranking candidates by score and off-target risk...", 84)
        results = self._rank_and_deduplicate(results, config)
        self._refine_top_results_with_rnaup(results, target_seq, transcriptome, config, report, 86, 93)
        self._annotate_decision_artifacts(results, target_seq_id, target_seq, transcriptome, config)
        results.sort(key=lambda x: x.get("recommendation_score", x["consensus_score"]), reverse=True)
        for i, r in enumerate(results):
            r["rank"] = i + 1

        report("Saving results...", 94)
        self._save_results(task_id, results, config)
        report("Design computation complete", 98)
        return results

    def _run_sgrna_mode(
        self,
        task_id: int,
        target_seq_id: str,
        target_seq: str,
        config: DesignConfig,
        transcriptome: TranscriptomeIndex,
        report: Callable,
    ) -> List[Dict]:
        """SpCas9 sgRNA mode: 20 nt spacer next to NGG PAM."""
        report("Scanning SpCas9 NGG sgRNA candidates...", 0)
        candidates = scan_sgrna_candidates(target_seq)
        report(f"Generated {len(candidates)} sgRNA candidates", 15)
        results = []
        for i, cand in enumerate(candidates):
            off_target = score_sgrna_offtargets(
                cand,
                transcriptome.sequences,
                exclude_target_id=target_seq_id,
            )
            cloning = design_sgrna_cloning_oligos(cand["spacer_dna"])
            genotyping = design_genotyping_primers(target_seq, cand["cut_site"])
            consensus_score = cand["on_target_score"]
            result_record = {
                "task_id": task_id,
                "rank": i + 1,
                "sequence": cand["guide_rna"],
                "position": f"{cand['position_start']}-{cand['position_end']}",
                "position_start": cand["position_start"],
                "position_end": cand["position_end"],
                "consensus_score": consensus_score,
                "passed": consensus_score >= 45 and off_target["passed"],
                "off_target": off_target,
                "sgrna": {
                    **cand,
                    "cloning_oligos": cloning,
                    "genotyping_primers": genotyping,
                },
                "primers": {
                    "sgrna_cloning_oligos": cloning,
                    "genotyping_primers": genotyping,
                },
            }
            results.append(result_record)
            if candidates:
                percent = 15 + (i + 1) / len(candidates) * 65
                report(f"Scoring sgRNA {i + 1}/{len(candidates)}", percent)

        report("Ranking sgRNAs by activity and off-target risk...", 84)
        results = self._rank_and_deduplicate(results, config)
        self._annotate_decision_artifacts(results, target_seq_id, target_seq, transcriptome, config)
        results.sort(key=lambda x: x.get("recommendation_score", x["consensus_score"]), reverse=True)
        for i, r in enumerate(results):
            r["rank"] = i + 1
        report("Saving results...", 94)
        self._save_results(task_id, results, config)
        report("Design computation complete", 98)
        return results

    def _annotate_recommendation_scores(self, results: List[Dict]):
        """Combine efficacy score and off-target risk into one sorting score."""
        for result in results:
            off_target = result.get("off_target") or {}
            risk_score = float(off_target.get("risk_score", 0) or 0)
            result["recommendation_score"] = round(float(result.get("consensus_score", 0) or 0) - risk_score * 0.45, 2)

    def _best_refinement_sequence(self, result: Dict) -> str:
        """Use the best Dicer product for pool-based designs, otherwise the candidate itself."""
        products = (result.get("pool") or {}).get("product_details", [])
        if products:
            best = max(products, key=lambda item: item.get("consensus_score", 0))
            return best.get("sequence", result.get("sequence", ""))
        return result.get("sequence", "")

    def _refine_top_results_with_rnaup(
        self,
        results: List[Dict],
        target_seq: str,
        transcriptome: TranscriptomeIndex,
        config: DesignConfig,
        report: Callable,
        start_percent: float,
        end_percent: float,
    ):
        """Run RNAup/explicit fallback only on top candidates and top risk transcripts."""
        top_n = int((config.thermodynamics or {}).get("rnaup_top_n", 0) or 0)
        if top_n <= 0 or not results:
            return
        ranked = sorted(results, key=lambda x: x.get("recommendation_score", x.get("consensus_score", 0)), reverse=True)
        selected = ranked[: min(top_n, len(ranked))]
        report(f"Running RNAup refinement for top {len(selected)} candidates...", start_percent)
        span = max(0.0, end_percent - start_percent)

        for i, result in enumerate(selected):
            query = self._best_refinement_sequence(result)
            risk_targets = (result.get("off_target") or {}).get("top_targets", [])
            risk_target_id = risk_targets[0]["target_id"] if risk_targets else None
            subject_seq = transcriptome.get_sequence(risk_target_id) if risk_target_id else target_seq
            subject_label = risk_target_id or "on_target"
            rnaup = self.thermo.rnaup(query, subject_seq) if query and subject_seq else None
            if rnaup:
                result["rnaup"] = {
                    **rnaup,
                    "target_id": subject_label,
                    "query_sequence": query,
                }
            if selected:
                report(
                    f"RNAup refinement {i + 1}/{len(selected)}",
                    start_percent + ((i + 1) / len(selected)) * span,
                )

    def _annotate_decision_artifacts(
        self,
        results: List[Dict],
        target_seq_id: str,
        target_seq: str,
        transcriptome: TranscriptomeIndex,
        config: DesignConfig,
    ):
        """Attach experiment-facing explanations, validation hits, maps and primers."""
        target_length = len(target_seq or "")
        for result in results:
            result["explanation"] = explain_result(result, mode=config.mode)
            result["validation_hits"] = build_validation_hits(result, transcriptome, max_hits=5)
            risk_positions = [
                hit.get("target_position", 0)
                for hit in result.get("validation_hits", [])[:3]
                if isinstance(hit.get("target_position", 0), int)
            ]
            result["region_map"] = render_region_map(
                target_length=target_length,
                start=result.get("position_start", 0),
                end=result.get("position_end", 0),
                risk_positions=risk_positions,
            )
            if config.mode == "long_dsRNA" and len(result.get("sequence", "")) >= 40:
                try:
                    primers = design_t7_primers(
                        result.get("sequence", ""),
                        result_id=f"candidate_{result.get('rank', '')}",
                    )
                    # Check primer specificity against transcriptome (exclude target itself)
                    if transcriptome and transcriptome.sequences:
                        from dsforge.core.primers import check_primer_specificity
                        exclude_ids = {target_seq_id}
                        fwd_spec = check_primer_specificity(
                            primers["forward_primer"]["sequence"],
                            transcriptome.sequences,
                            exclude_ids=exclude_ids,
                        )
                        rev_spec = check_primer_specificity(
                            primers["reverse_primer"]["sequence"],
                            transcriptome.sequences,
                            exclude_ids=exclude_ids,
                        )
                        primers["forward_primer"]["specificity"] = fwd_spec
                        primers["reverse_primer"]["specificity"] = rev_spec
                        if not fwd_spec["specific"] or not rev_spec["specific"]:
                            warning = (
                                "WARNING: Potential non-specific amplification detected. "
                                "Verify amplicon uniqueness with BLAST before ordering."
                            )
                            primers["notes"] = f"{primers.get('notes', '').strip()} {warning}".strip()
                    result["primers"] = primers
                except ValueError:
                    result["primers"] = {}

    def _save_results(self, task_id: int, results: List[Dict], config: DesignConfig):
        """保存结果到数据库"""
        with self.db._get_conn() as conn:
            for r in results:
                off_target = r.get("off_target") or {}
                top_targets = "; ".join(
                    target.get("target_id", "")
                    for target in off_target.get("top_targets", [])[:5]
                )
                cursor = conn.execute(
                    """
                    INSERT INTO results (
                        task_id, rank, candidate_seq, position_start, position_end,
                        consensus_score, passed_filters, risk_level, risk_score,
                        top_risk_targets, validation_direction, recommendation_score,
                        cluster_id, cluster_size, alternative_count, cluster_span,
                        explanation_json, validation_hits_json, primers_json, rnaup_json,
                        sgrna_json, region_map
                    )
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        task_id,
                        r["rank"],
                        r["sequence"],
                        r["position_start"],
                        r["position_end"],
                        r.get("consensus_score", 0),
                        1 if r.get("passed", False) else 0,
                        off_target.get("risk_level", "low"),
                        off_target.get("risk_score", 0),
                        top_targets,
                        off_target.get("validation_direction", ""),
                        r.get("recommendation_score", r.get("consensus_score", 0)),
                        r.get("cluster_id", r["rank"]),
                        r.get("cluster_size", 1),
                        r.get("alternative_count", 0),
                        r.get("cluster_span", r.get("position", "")),
                        json.dumps(r.get("explanation", {}), ensure_ascii=False),
                        json.dumps(r.get("validation_hits", []), ensure_ascii=False),
                        json.dumps(r.get("primers", {}), ensure_ascii=False),
                        json.dumps(r.get("rnaup", {}), ensure_ascii=False),
                        json.dumps(r.get("sgrna", {}), ensure_ascii=False),
                        r.get("region_map", ""),
                    ),
                )
                result_id = cursor.lastrowid

                for rule_name, rule_result in (r.get("rules") or {}).items():
                    conn.execute(
                        """
                        INSERT INTO rule_scores (result_id, rule_name, score, passed, violations)
                        VALUES (?, ?, ?, ?, ?)
                        """,
                        (
                            result_id,
                            rule_name,
                            rule_result.get("score", 0),
                            int(rule_result.get("passed", False)),
                            json.dumps(rule_result.get("violations", [])),
                        ),
                    )

                thermo = r.get("thermodynamics")
                if thermo:
                    conn.execute(
                        """
                        INSERT INTO thermodynamics (result_id, on_target_dg, seed_matches_count, high_risk_off_targets, rnaup_dg)
                        VALUES (?, ?, ?, ?, ?)
                        """,
                        (result_id, thermo.get("dg"), 0, 0, None),
                    )
                rnaup = r.get("rnaup")
                if rnaup:
                    conn.execute(
                        """
                        INSERT INTO thermodynamics (result_id, on_target_dg, seed_matches_count, high_risk_off_targets, rnaup_dg)
                        VALUES (?, ?, ?, ?, ?)
                        """,
                        (result_id, None, 0, 0, rnaup.get("dg")),
                    )

                for product in (r.get("pool") or {}).get("product_details", []):
                    conn.execute(
                        """
                        INSERT INTO pool_details (result_id, dicer_product_seq, cut_position, product_score)
                        VALUES (?, ?, ?, ?)
                        """,
                        (
                            result_id,
                            product.get("sequence", ""),
                            product.get("position"),
                            product.get("consensus_score"),
                        ),
                    )
            conn.commit()
