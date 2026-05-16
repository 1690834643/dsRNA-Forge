"""
Cached k-mer index for off-target risk ranking.
"""

import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from dsforge.core.sequence import DEFAULT_CACHE_DIR, normalize_sequence


VALID_BASES = set("AUGC")


def _is_clean_kmer(kmer: str) -> bool:
    return len(kmer) > 0 and set(kmer) <= VALID_BASES


class OffTargetRiskIndex:
    """7/16/20-mer transcriptome index for ranked off-target summaries."""

    def __init__(self, indexes: Dict[int, Dict], source_hash: Optional[str] = None):
        self.indexes = indexes
        self.source_hash = source_hash

    @classmethod
    def from_transcriptome(cls, transcriptome, cache_dir: Optional[Path] = None) -> "OffTargetRiskIndex":
        cache_dir = Path(cache_dir) if cache_dir is not None else getattr(transcriptome, "cache_dir", DEFAULT_CACHE_DIR)
        source_hash = getattr(transcriptome, "source_hash", None)
        cache_path = cache_dir / f"offtarget_{source_hash}_v3_7_16_20_27.json" if source_hash else None

        if cache_path and cache_path.exists():
            try:
                with open(cache_path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                if data.get("schema_version") != 3:
                    raise ValueError("Unsupported off-target cache schema")
                indexes = {int(k): v for k, v in data["indexes"].items()}
                return cls(indexes, source_hash=data.get("source_hash", source_hash))
            except Exception:
                pass

        index = cls.from_sequences(transcriptome.sequences, source_hash=source_hash)
        if cache_path:
            cache_dir.mkdir(parents=True, exist_ok=True)
            with open(cache_path, "w", encoding="utf-8") as f:
                json.dump(
                    {"indexes": index.indexes, "source_hash": source_hash, "schema_version": 3},
                    f,
                    ensure_ascii=False,
                )
        return index

    @classmethod
    def from_sequences(cls, sequences: Dict[str, str], source_hash: Optional[str] = None) -> "OffTargetRiskIndex":
        indexes = {7: defaultdict(dict), 16: defaultdict(dict), 20: defaultdict(dict), 27: defaultdict(dict)}
        for seq_id, raw_seq in sequences.items():
            seq = normalize_sequence(raw_seq)
            for k in indexes:
                if len(seq) < k:
                    continue
                for pos in range(len(seq) - k + 1):
                    kmer = seq[pos : pos + k]
                    if not _is_clean_kmer(kmer):
                        continue
                    hit = indexes[k][kmer].setdefault(seq_id, {"count": 0, "first_pos": pos})
                    hit["count"] += 1
        return cls({k: dict(v) for k, v in indexes.items()}, source_hash=source_hash)

    def assess_sequence(
        self,
        sequence: str,
        exclude_ids: Optional[Iterable[str]] = None,
        max_targets: int = 5,
        level_1_16bp: bool = True,
        level_2_20bp: bool = True,
        level_3_27bp: bool = False,
        seed_7nt: bool = True,
    ) -> Dict:
        """Return ranked risk details for one siRNA-like sequence."""
        seq = normalize_sequence(sequence)
        excluded = set(exclude_ids or [])
        per_target: Dict[str, Dict] = {}
        summary = {
            "level_1_16bp_hits": 0,
            "level_2_20bp_hits": 0,
            "level_3_27bp_hits": 0,
            "seed_7nt_hits": 0,
        }
        matches = []

        def add_hit(target_id: str, score: int, reason: str, match_type: str, length: int, position: int):
            if target_id in excluded:
                return
            target = per_target.setdefault(
                target_id,
                {"target_id": target_id, "risk_score": 0, "reasons": set(), "positions": []},
            )
            target["risk_score"] = max(target["risk_score"], score)
            target["reasons"].add(reason)
            target["positions"].append(position)
            matches.append({"target_id": target_id, "match_type": match_type, "length": length, "position": position})

        scan_levels = []
        if level_3_27bp:
            scan_levels.append((27, 100, "27bp near-perfect match", "27bp_consecutive", "level_3_27bp_hits"))
        if level_2_20bp:
            scan_levels.append((20, 95, "20bp continuous match", "20bp_consecutive", "level_2_20bp_hits"))
        if level_1_16bp:
            scan_levels.append((16, 75, "16bp continuous match", "16bp_consecutive", "level_1_16bp_hits"))

        for k, score, reason, match_type, summary_key in scan_levels:
            seen_targets = set()
            for pos in range(max(0, len(seq) - k + 1)):
                kmer = seq[pos : pos + k]
                for target_id, hit in self.indexes.get(k, {}).get(kmer, {}).items():
                    if target_id in excluded or target_id in seen_targets:
                        continue
                    seen_targets.add(target_id)
                    add_hit(target_id, score, reason, match_type, k, hit["first_pos"])
            summary[summary_key] = len(seen_targets)

        if seed_7nt:
            seed = seq[1:8] if len(seq) >= 8 else seq
            seed_targets = set()
            if _is_clean_kmer(seed):
                for target_id, hit in self.indexes.get(7, {}).get(seed, {}).items():
                    if target_id in excluded:
                        continue
                    seed_targets.add(target_id)
                    seed_score = min(55, 20 + hit["count"] * 2)
                    add_hit(target_id, seed_score, "7nt seed match", "seed_7nt", 7, hit["first_pos"])
            summary["seed_7nt_hits"] = len(seed_targets)

        top_targets = []
        for target in per_target.values():
            top_targets.append({
                "target_id": target["target_id"],
                "risk_score": target["risk_score"],
                "reasons": sorted(target["reasons"]),
                "first_position": min(target["positions"]) if target["positions"] else None,
            })
        top_targets.sort(key=lambda item: item["risk_score"], reverse=True)
        top_targets = top_targets[:max_targets]

        risk_score = top_targets[0]["risk_score"] if top_targets else 0
        if risk_score >= 80:
            risk_level = "high"
            validation = "优先避开该候选；如必须使用，先对 Top 风险转录本做 BLAST/qPCR 或片段比对验证。"
        elif risk_score >= 35 or summary["seed_7nt_hits"] > 5:
            risk_level = "medium"
            validation = "建议优先验证 Top 风险转录本，必要时换候选或缩小目标区域。"
        else:
            risk_level = "low"
            validation = "未发现强连续匹配；可作为低优先级脱靶验证对象。"

        return {
            "passed": risk_level != "high",
            "risk_level": risk_level,
            "risk_score": risk_score,
            "matches": matches,
            "summary": summary,
            "top_targets": top_targets,
            "risk_reasons": top_targets[0]["reasons"] if top_targets else [],
            "validation_direction": validation,
        }

    def assess_pool(
        self,
        products: List[Dict],
        exclude_ids: Optional[Iterable[str]] = None,
        level_1_16bp: bool = True,
        level_2_20bp: bool = True,
        level_3_27bp: bool = False,
        seed_7nt: bool = True,
        **_ignored,
    ) -> Dict:
        """Aggregate risk across Dicer products for long dsRNA/DsiRNA."""
        product_results = []
        high_risk = 0
        total_matches = 0
        target_scores: Dict[str, Dict] = {}

        for product in products:
            result = self.assess_sequence(
                product.get("sequence", ""),
                exclude_ids=exclude_ids,
                level_1_16bp=level_1_16bp,
                level_2_20bp=level_2_20bp,
                level_3_27bp=level_3_27bp,
                seed_7nt=seed_7nt,
            )
            if result["risk_level"] == "high":
                high_risk += 1
            total_matches += len(result["matches"])
            for target in result["top_targets"]:
                current = target_scores.get(target["target_id"])
                if current is None or target["risk_score"] > current["risk_score"]:
                    target_scores[target["target_id"]] = dict(target)
            product_results.append({**product, "off_target": result})

        top_targets = sorted(target_scores.values(), key=lambda item: item["risk_score"], reverse=True)[:5]
        max_score = top_targets[0]["risk_score"] if top_targets else 0
        risk_ratio = high_risk / len(products) if products else 0
        risk_score = min(100, round(max_score + risk_ratio * 10, 2))

        if risk_ratio > 0.2 or risk_score >= 80:
            pool_risk = "high"
            validation = "该 dsRNA pool 含高风险 Dicer 产物；优先验证 Top 风险转录本或换区段。"
        elif risk_ratio > 0.05 or risk_score >= 35:
            pool_risk = "medium"
            validation = "建议验证 Top 风险转录本，并优先选择风险分更低的区段。"
        else:
            pool_risk = "low"
            validation = "未发现强连续匹配；保留常规脱靶验证。"

        return {
            "pool_risk": pool_risk,
            "risk_level": pool_risk,
            "risk_score": risk_score,
            "passed": pool_risk != "high",
            "high_risk_products": high_risk,
            "total_matches": total_matches,
            "product_screening": product_results,
            "top_targets": top_targets,
            "risk_reasons": top_targets[0]["reasons"] if top_targets else [],
            "validation_direction": validation,
        }
