"""
脱靶筛查器

四级筛查层级：
Level 1: 16bp 连续匹配筛查 —— 排除高脱靶风险 dsRNA
Level 2: 20bp 连续匹配筛查 —— NTO 安全评估
Level 3: ViennaRNA RNAduplex 热力学评估 —— 种子区结合能风险预筛
Level 4: 种子区 7nt 匹配筛查 —— miRNA-like 脱靶
"""

from typing import List, Dict, Tuple, Optional, Iterable

from dsforge.core.offtarget_index import OffTargetRiskIndex


class OffTargetScreener:
    """脱靶筛查器"""

    def __init__(self, transcriptome_index):
        """
        Args:
            transcriptome_index: TranscriptomeIndex 实例
        """
        self.index = transcriptome_index
        self._thermo = None
        self._risk_index = None

    def prepare_index(self):
        """Build or load the cached off-target k-mer index."""
        if self._risk_index is None:
            self._risk_index = OffTargetRiskIndex.from_transcriptome(self.index)
        return self._risk_index

    def screen_sequence(
        self,
        sequence: str,
        exclude_ids: Optional[Iterable[str]] = None,
        level_1_16bp: bool = True,
        level_2_20bp: bool = True,
        level_3_27bp: bool = False,
        seed_7nt: bool = True,
        use_vienna: bool = True,
        seed_dg_threshold: float = -7.0,  # kcal/mol
        max_transcripts_for_thermo: int = 5000,
    ) -> Dict:
        """
        对单条 siRNA 序列进行脱靶筛查

        四级筛查 + ViennaRNA RNAduplex 热力学种子区评估：
        - Level 1: 16bp 连续匹配筛查
        - Level 2: 20bp 连续匹配筛查
        - Level 3: 27bp 近完美匹配
        - Level 4: 7nt 种子区匹配
        - ViennaRNA: RNAduplex 种子区热力学评估（如可用）

        Returns:
            {
                'passed': bool,
                'risk_level': str,      # 'low', 'medium', 'high'
                'matches': List[Dict],  # 匹配详情
                'summary': Dict,        # 各级统计
                'thermodynamics': Dict, # 热力学评估结果
            }
        """
        seq = sequence.upper().replace("T", "U")
        excluded = set(exclude_ids or [])
        self.prepare_index()
        risk = self._risk_index.assess_sequence(
            seq,
            exclude_ids=excluded,
            level_1_16bp=level_1_16bp,
            level_2_20bp=level_2_20bp,
            level_3_27bp=level_3_27bp,
            seed_7nt=seed_7nt,
        )

        thermo_results = {
            "vienna_available": False,
            "seed_hits": [],
            "min_dg": None,
        }

        # 获取种子区 (nt 2-8)
        seed = seq[1:8] if len(seq) >= 8 else seq

        # === ViennaRNA 热力学种子区评估 ===
        # 策略：先从缓存 k-mer 索引取出全部种子区/1nt 近似种子候选位点，
        # 再只对这些位点做 RNAduplex。这样不会按 FASTA 顺序截断整条转录本。
        if use_vienna:
            try:
                if self._thermo is None:
                    from dsforge.core.thermodynamics import ThermodynamicsCalculator
                    self._thermo = ThermodynamicsCalculator()
                calc = self._thermo
                if calc.available and len(seed) >= 6:
                    thermo_results["vienna_available"] = True
                    min_dg = 0.0
                    found = False

                    seed_sites = self._seed_site_candidates(seed, excluded)
                    thermo_results["candidate_seed_sites"] = len(seed_sites)
                    site_limit = int(max_transcripts_for_thermo or 0)
                    if site_limit > 0 and len(seed_sites) > site_limit:
                        thermo_results["truncated"] = True
                        thermo_results["unchecked_seed_hit_sites"] = len(seed_sites) - site_limit
                        thermo_results["truncation_reason"] = (
                            f"Checked {site_limit} indexed seed-hit sites. "
                            f"{len(seed_sites) - site_limit} additional seed-hit sites were not thermodynamically scored."
                        )
                        seed_sites = seed_sites[:site_limit]

                    # Step 2: 对可疑位置做 RNAduplex 热力学计算
                    for seq_id, i, target_sub, seed_mismatches in seed_sites:
                        result = calc.rnaduplex(seed, target_sub)
                        if result is not None:
                            dg = result["dg"]
                            if not found or dg < min_dg:
                                min_dg = dg
                                found = True
                            # 如果能量低于阈值，记录为热力学脱靶
                            if dg < seed_dg_threshold:
                                thermo_results["seed_hits"].append({
                                    "target_id": seq_id,
                                    "dg": round(dg, 2),
                                    "position": i,
                                    "seed": seed,
                                    "target_seq": target_sub,
                                    "seed_mismatches": seed_mismatches,
                                    "structure": result.get("structure", ""),
                                })

                    if found:
                        thermo_results["min_dg"] = round(min_dg, 2)
            except Exception as e:
                # ViennaRNA 不可用或出错，降级为纯序列比对
                thermo_results["error"] = str(e)

        return {
            "passed": risk["passed"],
            "risk_level": risk["risk_level"],
            "risk_score": risk["risk_score"],
            "matches": risk["matches"],
            "summary": risk["summary"],
            "thermodynamics": thermo_results,
            "top_targets": risk["top_targets"],
            "risk_reasons": risk["risk_reasons"],
            "validation_direction": risk["validation_direction"],
        }

    def _seed_site_candidates(self, seed: str, excluded: set) -> List[Tuple[str, int, str, int]]:
        """Return indexed seed and 1-mismatch seed sites in deterministic risk order."""
        seed = seed.upper().replace("T", "U")
        if len(seed) < 6 or set(seed) - set("AUGC"):
            return []

        k_index = self._risk_index.indexes.get(len(seed), {})
        variants = self._seed_variants(seed)
        candidates: Dict[Tuple[str, int, str], int] = {}
        for variant, mismatches in variants:
            for seq_id in k_index.get(variant, {}):
                if seq_id in excluded:
                    continue
                target_seq = self.index.sequences.get(seq_id, "").upper().replace("T", "U")
                pos = target_seq.find(variant)
                while pos != -1:
                    key = (seq_id, pos, variant)
                    current = candidates.get(key)
                    if current is None or mismatches < current:
                        candidates[key] = mismatches
                    pos = target_seq.find(variant, pos + 1)

        rows = [
            (seq_id, pos, target_sub, mismatches)
            for (seq_id, pos, target_sub), mismatches in candidates.items()
        ]
        rows.sort(key=lambda item: (item[3], item[0], item[1], item[2]))
        return rows

    def _seed_variants(self, seed: str) -> List[Tuple[str, int]]:
        """Return exact seed and all one-mismatch seed variants."""
        variants = [(seed, 0)]
        bases = "AUGC"
        for pos, base in enumerate(seed):
            for replacement in bases:
                if replacement == base:
                    continue
                variants.append((seed[:pos] + replacement + seed[pos + 1 :], 1))
        return variants

    def _has_consecutive_match(self, query: str, target: str, min_len: int) -> bool:
        """检查 query 是否有长度 >= min_len 的连续子串匹配到 target"""
        if len(query) < min_len:
            return False
        for i in range(len(query) - min_len + 1):
            substring = query[i : i + min_len]
            if substring in target:
                return True
        return False

    def _has_near_perfect_match(
        self, query: str, target: str, min_len: int, max_mismatches: int = 2
    ) -> bool:
        """检查 query 是否有近完美匹配（允许少量错配）"""
        if len(query) < min_len:
            return False
        for i in range(len(query) - min_len + 1):
            substring = query[i : i + min_len]
            for j in range(len(target) - min_len + 1):
                target_sub = target[j : j + min_len]
                mismatches = sum(a != b for a, b in zip(substring, target_sub))
                if mismatches <= max_mismatches:
                    return True
        return False

    def screen_pool(
        self,
        products: List[Dict],
        exclude_ids: Optional[Iterable[str]] = None,
        **kwargs,
    ) -> Dict:
        """
        对 siRNA Pool 进行脱靶筛查

        Returns:
            {
                'pool_risk': str,        # 'low', 'medium', 'high'
                'high_risk_products': int,
                'total_matches': int,
                'product_screening': List[Dict],
            }
        """
        self.prepare_index()
        return self._risk_index.assess_pool(products, exclude_ids=exclude_ids, **kwargs)
