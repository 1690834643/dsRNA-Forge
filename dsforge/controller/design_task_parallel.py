"""
多进程设计任务控制器
集成 ParallelScheduler 到设计任务中
"""

from typing import Dict, List, Callable
import multiprocessing as mp
from dataclasses import asdict

from dsforge.controller.design_task import DesignCancelled, DesignTask, DesignConfig
from dsforge.core.sequence import TranscriptomeIndex
from dsforge.core.scoring.base import evaluate_all_rules
from dsforge.core.dicer import predict_dicer_products, calculate_pool_score
from dsforge.core.offtarget import OffTargetScreener
from dsforge.core.offtarget_index import OffTargetRiskIndex
from dsforge.core.thermodynamics import ThermodynamicsCalculator
from dsforge.database.manager import DatabaseManager


_WORKER_RISK_INDEX = None


def _init_worker_risk_index(risk_index):
    """Install a pickled risk index once per worker process."""
    global _WORKER_RISK_INDEX
    _WORKER_RISK_INDEX = risk_index


def _score_sirna_batch(args: dict):
    """
    子进程函数：评分一批 siRNA 候选
    必须作为模块级函数以便 pickle
    """
    batch = args["batch"]
    target_seq = args["target_seq"]
    target_seq_id = args["target_seq_id"]
    enabled_rules = args["enabled_rules"]
    off_target_config = args["off_target_config"]
    transcriptome_seqs = args.get("transcriptome_seqs")

    if transcriptome_seqs is not None:
        risk_index = OffTargetRiskIndex.from_sequences(transcriptome_seqs)
    else:
        risk_index = _WORKER_RISK_INDEX

    results = []
    for cand in batch:
        seq = cand["sequence"]
        rule_results = evaluate_all_rules(seq, enabled_rules)
        consensus = rule_results.get("consensus", {})

        off_target = {"passed": True, "risk_level": "low", "risk_score": 0, "matches": [], "top_targets": []}
        if risk_index is not None:
            off_target = risk_index.assess_sequence(
                seq,
                exclude_ids={target_seq_id},
                **off_target_config,
            )

        results.append({
            "sequence": seq,
            "start": cand["start"],
            "end": cand["end"],
            "consensus_score": consensus.get("score", 0),
            "passed": consensus.get("passed", False) and off_target["passed"],
            "rules": rule_results,
            "off_target": off_target,
        })
    return results


def _score_long_dsRNA_batch(args: dict):
    """
    子进程函数：评分一批长 dsRNA 区段
    """
    batch = args["batch"]
    enabled_rules = args["enabled_rules"]
    off_target_config = args["off_target_config"]
    target_seq_id = args["target_seq_id"]
    transcriptome_seqs = args.get("transcriptome_seqs")
    risk_index_override = args.get("risk_index")

    if transcriptome_seqs is not None:
        risk_index = risk_index_override or OffTargetRiskIndex.from_sequences(transcriptome_seqs)
    else:
        risk_index = _WORKER_RISK_INDEX

    results = []
    for cand in batch:
        seq = cand["sequence"]
        products = predict_dicer_products(seq)
        pool_result = calculate_pool_score(
            products, evaluate_all_rules, enabled_rules
        )
        off_target = {"passed": True, "risk_level": "low", "risk_score": 0, "top_targets": []}
        if risk_index is not None:
            off_target = risk_index.assess_pool(
                products,
                exclude_ids={target_seq_id},
                **off_target_config,
            )

        results.append({
            "sequence": seq,
            "start": cand["start"],
            "end": cand["end"],
            "consensus_score": pool_result["pool_score"],
            "passed": pool_result["pool_score"] > 50 and off_target["passed"],
            "pool": pool_result,
            "off_target": off_target,
        })
    return results


class ParallelDesignTask(DesignTask):
    """
    多进程设计任务控制器
    继承单线程版本，添加多进程支持
    """

    def __init__(self, db_manager=None, n_cores: int = -1, batch_size: int = 100):
        super().__init__(db_manager)
        self.n_cores = n_cores if n_cores > 0 else max(1, mp.cpu_count() - 1)
        self.batch_size = batch_size

    def _run_sirna_mode_parallel(
        self,
        task_id: int,
        target_seq_id: str,
        target_seq: str,
        config: DesignConfig,
        transcriptome: TranscriptomeIndex,
        report: Callable,
    ) -> List[Dict]:
        """siRNA 模式 — 多进程评分"""
        from dsforge.core.sequence import generate_candidates

        report("Generating siRNA candidates...", 0)
        self._ensure_candidate_budget(target_seq, 21, 21, config)
        candidates = list(generate_candidates(
            target_seq, "siRNA", 21, 21,
            gc_min=config.gc_min, gc_max=config.gc_max,
            exclude_poly=config.exclude_poly_n,
            max_candidates=config.max_raw_candidates,
        ))
        report(f"Generated {len(candidates)} candidates", 10)

        # 分批
        batches = [
            candidates[i:i + self.batch_size]
            for i in range(0, len(candidates), self.batch_size)
        ]

        # 准备参数（picklable）
        off_target_config = config.off_target_levels or {}
        report("Building or reusing off-target risk index...", 12)
        risk_index = OffTargetRiskIndex.from_transcriptome(transcriptome) if transcriptome else None

        args = [
            {
                "batch": batch,
                "target_seq": target_seq,
                "target_seq_id": target_seq_id,
                "enabled_rules": config.enabled_rules,
                "off_target_config": off_target_config,
            }
            for batch in batches
        ]

        report(f"Scoring with {self.n_cores} cores...", 15)
        results = []

        with mp.Pool(processes=self.n_cores, initializer=_init_worker_risk_index, initargs=(risk_index,)) as pool:
            for i, batch_results in enumerate(pool.imap_unordered(_score_sirna_batch, args)):
                results.extend(batch_results)
                percent = 15 + (i + 1) / len(batches) * 60
                report(f"Batch {i+1}/{len(batches)} done", percent)

        report("Ranking candidates by score and off-target risk...", 82)
        full_results = []
        for i, r in enumerate(results):
            full_results.append({
                "task_id": task_id,
                "rank": i + 1,
                "sequence": r["sequence"],
                "position": f"{r['start']}-{r['end']}",
                "position_start": r["start"],
                "position_end": r["end"],
                "consensus_score": r["consensus_score"],
                "passed": r["passed"],
                "rules": r["rules"],
                "off_target": r["off_target"],
            })
        full_results = self._rank_and_deduplicate(full_results, config)

        self._refine_top_results_with_rnaup(full_results, target_seq, transcriptome, config, report, 86, 93)
        self._annotate_decision_artifacts(full_results, target_seq_id, target_seq, transcriptome, config)
        full_results.sort(key=lambda x: x.get("recommendation_score", x["consensus_score"]), reverse=True)
        for i, r in enumerate(full_results):
            r["rank"] = i + 1

        report("Saving results...", 94)
        self._save_results(task_id, full_results, config)
        report("Design computation complete", 98)
        return full_results

    def _run_long_dsRNA_mode_parallel(
        self,
        task_id: int,
        target_seq_id: str,
        target_seq: str,
        config: DesignConfig,
        transcriptome: TranscriptomeIndex,
        report: Callable,
    ) -> List[Dict]:
        """长 dsRNA 模式 — 多进程区段评估"""
        from dsforge.core.sequence import generate_candidates

        report("Generating long dsRNA candidates...", 0)
        self._ensure_candidate_budget(target_seq, config.length_min, config.length_max, config)
        candidates = list(generate_candidates(
            target_seq, "long_dsRNA",
            min_len=config.length_min, max_len=config.length_max,
            gc_min=config.gc_min, gc_max=config.gc_max,
            exclude_poly=config.exclude_poly_n,
            max_candidates=config.max_raw_candidates,
        ))
        report(f"Generated {len(candidates)} candidate regions", 5)

        batches = [
            candidates[i:i + self.batch_size]
            for i in range(0, len(candidates), self.batch_size)
        ]

        off_target_config = config.off_target_levels or {}
        report("Building or reusing off-target risk index...", 8)
        risk_index = OffTargetRiskIndex.from_transcriptome(transcriptome) if transcriptome else None
        args = [
            {
                "batch": batch,
                "enabled_rules": config.enabled_rules,
                "off_target_config": off_target_config,
                "target_seq_id": target_seq_id,
            }
            for batch in batches
        ]

        report(f"Evaluating with {self.n_cores} cores...", 10)
        results = []

        with mp.Pool(processes=self.n_cores, initializer=_init_worker_risk_index, initargs=(risk_index,)) as pool:
            for i, batch_results in enumerate(pool.imap_unordered(_score_long_dsRNA_batch, args)):
                results.extend(batch_results)
                percent = 10 + (i + 1) / len(batches) * 65
                report(f"Batch {i+1}/{len(batches)} done", percent)

        report("Ranking candidates by score and off-target risk...", 82)
        full_results = []
        for i, r in enumerate(results):
            full_results.append({
                "task_id": task_id,
                "rank": i + 1,
                "sequence": r["sequence"],
                "position": f"{r['start']}-{r['end']}",
                "position_start": r["start"],
                "position_end": r["end"],
                "consensus_score": r["consensus_score"],
                "passed": r["passed"],
                "pool": r["pool"],
                "off_target": r["off_target"],
            })
        full_results = self._rank_and_deduplicate(full_results, config)

        self._refine_top_results_with_rnaup(full_results, target_seq, transcriptome, config, report, 86, 93)
        self._annotate_decision_artifacts(full_results, target_seq_id, target_seq, transcriptome, config)
        full_results.sort(key=lambda x: x.get("recommendation_score", x["consensus_score"]), reverse=True)
        for i, r in enumerate(full_results):
            r["rank"] = i + 1

        report("Saving results...", 94)
        self._save_results(task_id, full_results, config)
        report("Design computation complete", 98)
        return full_results

    def run_parallel(
        self,
        transcriptome: TranscriptomeIndex,
        target_seq_id: str,
        config: DesignConfig,
        progress_callback: Callable = None,
    ) -> Dict:
        """
        多进程执行设计任务
        """
        target_seq = transcriptome.get_sequence(target_seq_id)
        if target_seq is None:
            raise ValueError(f"Target sequence '{target_seq_id}' not found")

        task_id = self.db.create_task(
            mode=config.mode,
            target_seq_id=target_seq_id,
            target_seq=target_seq,
            params={**asdict(config), "n_cores": self.n_cores, "batch_size": self.batch_size},
        )
        self.db.update_task_status(task_id, "running")

        def _report(step: str, percent: float):
            if progress_callback:
                progress_callback(step, percent)

        try:
            if config.mode == "siRNA":
                results = self._run_sirna_mode_parallel(task_id, target_seq_id, target_seq, config, transcriptome, _report)
            elif config.mode == "long_dsRNA":
                results = self._run_long_dsRNA_mode_parallel(task_id, target_seq_id, target_seq, config, transcriptome, _report)
            elif config.mode == "sgRNA":
                results = self._run_sgrna_mode(task_id, target_seq_id, target_seq, config, transcriptome, _report)
            else:
                # DsiRNA 和单线程模式回退
                results = self._run_dsiRNA_mode(task_id, target_seq_id, target_seq, config, transcriptome, _report)

            self.db.update_task_status(task_id, "completed")

            return {
                "task_id": task_id,
                "mode": config.mode,
                "target_seq": target_seq,
                "config": config,
                "results": results,
                "summary": {
                    "total_candidates": len(results),
                    "passed_candidates": sum(1 for r in results if r.get("passed", False)),
                    "raw_candidates": getattr(self, "_last_raw_candidate_count", len(results)),
                    "nonredundant_candidates": len(results),
                    "deduplicated_candidates": max(0, getattr(self, "_last_raw_candidate_count", len(results)) - len(results)),
                },
            }

        except DesignCancelled:
            self.db.update_task_status(task_id, "cancelled")
            raise
        except Exception as e:
            self.db.update_task_status(task_id, "failed")
            raise
