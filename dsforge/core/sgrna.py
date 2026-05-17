"""
SpCas9 sgRNA design helpers.

This module implements an offline, conservative subset of mainstream sgRNA
design workflows: 20 nt spacer discovery next to NGG PAMs, basic on-target
heuristics, mismatch-weighted off-target ranking, cloning oligos, and simple
genotyping PCR primer suggestions.
"""

from typing import Dict, Iterable, List, Optional

from dsforge.core.primers import design_t7_primers
from dsforge.core.sequence import gc_content, has_poly_repeat, normalize_sequence


DNA_COMPLEMENT = str.maketrans("ACGTU", "TGCAA")
VALID_DNA = set("ACGT")
SGRNA_MISMATCH_BUCKETS = tuple(f"{idx}M" for idx in range(6))
SGRNA_SEED_INDEX_CHUNKS = ((0, 4), (4, 8), (8, 11), (11, 14), (14, 17), (17, 20))
STOP_CODONS = {"TAA", "TAG", "TGA"}
RESTRICTION_ENZYME_SITES = {
    "BsaI": ("GGTCTC", "GAGACC"),
    "BbsI": ("GAAGAC", "GTCTTC"),
    "BsmBI": ("CGTCTC", "GAGACG"),
}


class SgRNAReferenceSites:
    """Precomputed Cas9 sites plus a no-false-negative chunk index."""

    def __init__(self, sites: List[Dict]):
        self.sites = sites
        self.chunk_index: Dict[tuple, List[int]] = {}
        for idx, site in enumerate(sites):
            spacer = site.get("spacer_dna", "")
            if len(spacer) != 20:
                continue
            for chunk_id, (start, end) in enumerate(SGRNA_SEED_INDEX_CHUNKS):
                self.chunk_index.setdefault((chunk_id, spacer[start:end]), []).append(idx)

    def __iter__(self):
        return iter(self.sites)

    def __len__(self):
        return len(self.sites)

    def __getitem__(self, idx):
        return self.sites[idx]

    def candidate_sites(self, spacer: str, max_mismatches: int) -> List[Dict]:
        if max_mismatches >= len(SGRNA_SEED_INDEX_CHUNKS):
            return self.sites
        site_indices = set()
        for chunk_id, (start, end) in enumerate(SGRNA_SEED_INDEX_CHUNKS):
            site_indices.update(self.chunk_index.get((chunk_id, spacer[start:end]), ()))
        return [self.sites[idx] for idx in site_indices]


def to_dna(sequence: str) -> str:
    return normalize_sequence(sequence).replace("U", "T")


def reverse_complement_dna(sequence: str) -> str:
    return to_dna(sequence).translate(DNA_COMPLEMENT)[::-1]


def _valid_spacer(spacer: str) -> bool:
    return len(spacer) == 20 and set(spacer) <= VALID_DNA


def _pam_matches(pam_seq: str, pattern: str) -> bool:
    pam_seq = to_dna(pam_seq)
    pattern = pattern.upper()
    if len(pam_seq) != 3:
        return False
    if pattern == "NGG":
        return pam_seq[1:] == "GG"
    if pattern == "NRG":
        return pam_seq[1] in {"A", "G"} and pam_seq[2] == "G"
    raise ValueError("Only SpCas9 NGG and NRG PAM patterns are supported")


def _scan_spcas9_sites(sequence: str, pam: str = "NGG") -> List[Dict]:
    """Scan both strands for guide-oriented SpCas9 sites matching a PAM pattern."""
    dna = to_dna(sequence)
    candidates: List[Dict] = []

    for i in range(20, len(dna) - 2):
        pam_seq = dna[i : i + 3]
        if _pam_matches(pam_seq, pam):
            spacer = dna[i - 20 : i]
            if _valid_spacer(spacer):
                score = _score_spacer(spacer)
                cut_site = i - 3
                candidates.append({
                    "spacer_dna": spacer,
                    "guide_rna": spacer.replace("T", "U"),
                    "pam": pam_seq,
                    "genomic_pam": pam_seq,
                    "strand": "+",
                    "position_start": i - 20,
                    "position_end": i,
                    "pam_start": i,
                    "pam_end": i + 3,
                    "cut_site": cut_site,
                    **score,
                })

    # Reverse-strand target appears as the reverse complement of the guide PAM
    # on the provided plus-strand sequence: CCN for NGG and CYN for NRG.
    for i in range(0, len(dna) - 22):
        genomic_pam = dna[i : i + 3]
        guide_oriented_pam = reverse_complement_dna(genomic_pam)
        if _pam_matches(guide_oriented_pam, pam):
            protospacer_on_plus = dna[i + 3 : i + 23]
            spacer = reverse_complement_dna(protospacer_on_plus)
            if _valid_spacer(spacer):
                score = _score_spacer(spacer)
                cut_site = i + 6
                candidates.append({
                    "spacer_dna": spacer,
                    "guide_rna": spacer.replace("T", "U"),
                    "pam": guide_oriented_pam,
                    "genomic_pam": genomic_pam,
                    "strand": "-",
                    "position_start": i + 3,
                    "position_end": i + 23,
                    "pam_start": i,
                    "pam_end": i + 3,
                    "cut_site": cut_site,
                    **score,
                })

    return candidates


def _score_spacer(spacer: str) -> Dict:
    gc = gc_content(spacer)
    score = 100.0
    notes = []
    if gc < 35:
        penalty = min(30, (35 - gc) * 1.4)
        score -= penalty
        notes.append("GC below preferred 35-65% range")
    elif gc > 65:
        penalty = min(30, (gc - 65) * 1.4)
        score -= penalty
        notes.append("GC above preferred 35-65% range")
    if has_poly_repeat(spacer.replace("T", "U"), 4)[0]:
        score -= 18
        notes.append("Contains homopolymer run")
    if "TTTT" in spacer:
        score -= 22
        notes.append("Contains U6 terminator-like TTTT motif")
    if spacer[0] == "G":
        score += 4
        notes.append("Starts with G, convenient for U6 expression")
    elif spacer[0] in {"A", "T"}:
        score -= 4
    # Avoid extreme base imbalance in PAM-proximal seed region.
    seed_gc = gc_content(spacer[-12:])
    if seed_gc < 25 or seed_gc > 80:
        score -= 8
        notes.append("PAM-proximal seed GC is extreme")
    return {
        "on_target_score": round(max(0.0, min(100.0, score)), 2),
        "gc_percent": round(gc, 1),
        "notes": notes,
    }


def scan_sgrna_candidates(sequence: str, pam: str = "NGG") -> List[Dict]:
    """Scan both strands for SpCas9 NGG candidates."""
    candidates = _scan_spcas9_sites(sequence, pam=pam)
    candidates.sort(key=lambda item: (item["on_target_score"], -abs(50 - item["gc_percent"])), reverse=True)
    for rank, candidate in enumerate(candidates, start=1):
        candidate["rank"] = rank
    return candidates


def _find_orfs(dna: str, min_orf_nt: int = 90) -> List[Dict]:
    """Return simple ATG-to-stop ORF candidates on the provided transcript strand."""
    orfs = []
    dna_len = len(dna)
    for frame in range(3):
        pos = frame
        while pos <= dna_len - 3:
            if dna[pos : pos + 3] != "ATG":
                pos += 3
                continue
            stop_pos = None
            scan = pos + 3
            while scan <= dna_len - 3:
                if dna[scan : scan + 3] in STOP_CODONS:
                    stop_pos = scan + 3
                    break
                scan += 3
            end = stop_pos if stop_pos is not None else dna_len - ((dna_len - pos) % 3)
            if end - pos >= min_orf_nt:
                orfs.append({
                    "start": pos,
                    "end": end,
                    "frame": frame,
                    "length": end - pos,
                    "has_stop": stop_pos is not None,
                })
            pos += 3
    return orfs


def prepare_sgrna_design_sequence(sequence: str, min_orf_nt: int = 90) -> Dict:
    """Choose the sequence used for sgRNA design and record CDS/mRNA advice.

    If a transcript-like input contains a clear ATG-to-stop ORF, guides are
    designed on that inferred CDS slice and coordinates are mapped back to the
    original input. If no confident ORF is found, the full input is retained.
    """
    dna = to_dna(sequence)
    advice = [
        "sgRNA 设计建议优先提交 CDS 序列；如果输入为 mRNA/cDNA，本工具会尝试用最长 ATG-to-stop ORF 推断 CDS。",
    ]
    result = {
        "design_sequence": dna,
        "source_start": 0,
        "source_end": len(dna),
        "source_length": len(dna),
        "input_type": "target_fragment",
        "cds_inferred": False,
        "uses_cds": False,
        "orf_frame": None,
        "orf_has_stop": False,
        "advice": advice,
    }
    if len(dna) < 23 or set(dna) - VALID_DNA:
        advice.append("未执行 CDS 推断：序列过短或含 N/模糊碱基；请尽量提供 A/C/G/T 明确的 CDS。")
        return result

    orfs = _find_orfs(dna, min_orf_nt=min_orf_nt)
    if not orfs:
        advice.append("未检测到可信 ATG-to-stop CDS；当前按完整输入片段设计，结果需要人工确认是否位于编码区。")
        return result

    best = max(orfs, key=lambda item: (item["length"], item["has_stop"]))
    result.update({
        "design_sequence": dna[best["start"] : best["end"]],
        "source_start": best["start"],
        "source_end": best["end"],
        "input_type": "cds" if best["start"] == 0 and best["end"] == len(dna) else "mrna_inferred_cds",
        "cds_inferred": best["start"] != 0 or best["end"] != len(dna),
        "uses_cds": True,
        "orf_frame": best["frame"],
        "orf_has_stop": best["has_stop"],
    })
    if result["cds_inferred"]:
        advice.append(
            f"已从输入 mRNA/cDNA 推断 CDS：原始坐标 {best['start']}-{best['end']}，sgRNA 坐标会映射回原始序列。"
        )
    else:
        advice.append("输入看起来像 CDS；排序会优先推荐 CDS 前段候选。")
    advice.append("Cas9 knockout 通常优先选择靠前 CDS 外显子/前段 CDS 位点，避开仅位于 UTR 的候选。")
    return result


def annotate_sgrna_cds_priority(candidates: List[Dict], design_info: Dict) -> List[Dict]:
    """Annotate candidates with CDS-aware coordinates and early-CDS priority."""
    design_length = max(1, len(design_info.get("design_sequence", "")) or 1)
    offset = int(design_info.get("source_start", 0) or 0)
    uses_cds = bool(design_info.get("uses_cds"))
    annotated = []
    for candidate in candidates:
        item = dict(candidate)
        cut_site = int(item.get("cut_site", item.get("position_start", 0)) or 0)
        fraction = max(0.0, min(1.0, cut_site / design_length))
        if fraction <= 0.35:
            region = "front_cds" if uses_cds else "front_input"
            bonus = 8.0 if uses_cds else 0.0
        elif fraction <= 0.65:
            region = "middle_cds" if uses_cds else "middle_input"
            bonus = 3.0 if uses_cds else 0.0
        else:
            region = "late_cds" if uses_cds else "late_input"
            bonus = -8.0 if uses_cds else 0.0
        item.update({
            "source_position_start": offset + int(item.get("position_start", 0) or 0),
            "source_position_end": offset + int(item.get("position_end", 0) or 0),
            "source_pam_start": offset + int(item.get("pam_start", item.get("position_end", 0)) or 0),
            "source_pam_end": offset + int(item.get("pam_end", item.get("position_end", 0)) or 0),
            "source_cut_site": offset + cut_site,
            "cds_position_percent": round(fraction * 100, 1),
            "cds_region": region,
            "cds_priority_bonus": bonus,
            "locus_priority_score": round(float(item.get("on_target_score", 0) or 0) + bonus, 2),
            "design_input_type": design_info.get("input_type", "target_fragment"),
            "cds_source_start": design_info.get("source_start", 0),
            "cds_source_end": design_info.get("source_end", design_info.get("source_length", 0)),
            "input_advice": list(design_info.get("advice", [])),
        })
        annotated.append(item)
    annotated.sort(
        key=lambda item: (
            item.get("locus_priority_score", item.get("on_target_score", 0)),
            item.get("on_target_score", 0),
            -abs(50 - item.get("gc_percent", 50)),
        ),
        reverse=True,
    )
    for rank, item in enumerate(annotated, start=1):
        item["rank"] = rank
    return annotated


def _mismatch_positions(a: str, b: str) -> List[int]:
    return [idx + 1 for idx, (x, y) in enumerate(zip(a, b)) if x != y]


def _cfd_like_score(mismatch_positions: Iterable[int]) -> float:
    score = 100.0
    for pos in mismatch_positions:
        # PAM-proximal seed positions receive stronger penalties.
        if pos >= 13:
            score *= 0.18
        elif pos >= 8:
            score *= 0.38
        else:
            score *= 0.62
    return round(score, 2)


def _empty_mismatch_counts() -> Dict[str, int]:
    return {bucket: 0 for bucket in SGRNA_MISMATCH_BUCKETS}


def _count_mismatch_buckets(hits: List[Dict], seed_only: bool = False) -> Dict[str, int]:
    counts = _empty_mismatch_counts()
    for hit in hits:
        if seed_only and not hit.get("seed_12_match"):
            continue
        mismatches = int(hit.get("mismatches", 0) or 0)
        if mismatches <= 5:
            counts[f"{mismatches}M"] += 1
    return counts


def _sgrnacas9_risk_label(summary: Dict) -> str:
    counts = summary["mismatch_counts"]
    pot_counts = summary["pot_mismatch_counts"]
    if counts["0M"] >= 2:
        return "Repeat_sites_or_bad?"
    if counts["0M"] == 1:
        return "High_risk"
    if counts["1M"] or pot_counts["1M"]:
        return "High_risk"
    if counts["2M"] or pot_counts["2M"]:
        return "Moderate_risk"
    if pot_counts["3M"] or pot_counts["4M"] or pot_counts["5M"]:
        return "Low_risk"
    return "Best"


def _risk_from_hits(hits: List[Dict], max_mismatches: int = 5) -> Dict:
    summary = {
        "mismatch_counts": _count_mismatch_buckets(hits),
        "pot_mismatch_counts": _count_mismatch_buckets(hits, seed_only=True),
        "total_ot": len(hits),
        "total_pot": sum(1 for h in hits if h.get("seed_12_match")),
        "nrg_pam_hits": len(hits),
        "nag_pam_hits": sum(1 for h in hits if h.get("pam", "")[1:2] == "A"),
        "alternative_pam_hits": sum(1 for h in hits if h.get("pam", "")[1:2] == "A"),
        "max_mismatches_checked": max_mismatches,
        "reference_scope": "current_reference_sequences",
        "reference_scope_note": "sgRNA off-target scan covers only loaded transcriptome/reference/background FASTA sequences; use genome FASTA for genome-scale Cas9 off-target screening.",
    }
    summary.update({
        "perfect_hits": summary["mismatch_counts"]["0M"],
        "one_mismatch_hits": summary["mismatch_counts"]["1M"],
        "two_mismatch_hits": summary["mismatch_counts"]["2M"],
        "three_mismatch_hits": summary["mismatch_counts"]["3M"],
        "four_mismatch_hits": summary["mismatch_counts"]["4M"],
        "five_mismatch_hits": summary["mismatch_counts"]["5M"],
    })
    summary["sgrnacas9_risk_evaluation"] = _sgrnacas9_risk_label(summary)

    if not hits:
        return {
            "passed": True,
            "risk_level": "low",
            "risk_score": 0,
            "top_targets": [],
            "matches": [],
            "summary": summary,
            "risk_reasons": [],
            "validation_direction": f"未发现 <={max_mismatches} mismatch 且邻近 SpCas9 NRG PAM 的强 off-target；保留常规验证。",
        }

    hits.sort(key=lambda item: (item["risk_score"], -item["mismatches"]), reverse=True)
    top = hits[:5]
    max_score = top[0]["risk_score"]
    seed_pot_count = summary["total_pot"]
    if any(h["mismatches"] == 0 for h in top) or max_score >= 85:
        level = "high"
        validation = "存在高风险 Cas9 off-target；优先换 sgRNA，或对 Top 位点做扩增测序/ICE/TIDE 验证。"
    elif max_score >= 35 or any(h["mismatches"] <= 1 for h in top):
        level = "medium"
        validation = "建议对 Top off-target 位点设计 PCR 引物并做靶向测序验证。"
    elif seed_pot_count:
        level = "medium"
        validation = "存在 12 nt seed 完全匹配的 POT 位点；建议优先对这些 Top 位点做 PCR 扩增测序验证。"
    else:
        level = "low"
        validation = "仅发现低分 off-target；保留常规靶点验证。"

    return {
        "passed": level != "high",
        "risk_level": level,
        "risk_score": max_score,
        "top_targets": top,
        "matches": [
            {
                "target_id": h["target_id"],
                "match_type": f"Cas9_{h['mismatches']}mm",
                "length": 20,
                "position": h["position"],
                "strand": h["strand"],
                "pam": h.get("pam", ""),
                "mismatch_positions": h.get("mismatch_positions", []),
                "seed_12_match": h.get("seed_12_match", False),
                "target_spacer": h.get("target_spacer", ""),
                "target_protospacer_pam": h.get("target_protospacer_pam", ""),
                "validation_window_start": h.get("validation_window_start", 0),
                "validation_window_end": h.get("validation_window_end", 0),
            }
            for h in top
        ] + [
            {
                "target_id": h["target_id"],
                "match_type": f"Cas9_seed12_POT_{h['mismatches']}mm",
                "length": 12,
                "position": h["position"],
                "strand": h["strand"],
                "pam": h.get("pam", ""),
                "mismatch_positions": h.get("mismatch_positions", []),
                "seed_12_match": True,
                "target_spacer": h.get("target_spacer", ""),
                "target_protospacer_pam": h.get("target_protospacer_pam", ""),
                "validation_window_start": h.get("validation_window_start", 0),
                "validation_window_end": h.get("validation_window_end", 0),
            }
            for h in top
            if h.get("seed_12_match")
        ],
        "summary": summary,
        "risk_reasons": top[0].get("reasons", []),
        "validation_direction": validation,
    }


def build_sgrna_reference_sites(
    reference_sequences: Dict[str, str],
    pam: str = "NRG",
    exclude_target_ids: Optional[Iterable[str]] = None,
) -> SgRNAReferenceSites:
    """Pre-scan reference sequences into guide-oriented SpCas9 sites."""
    excluded = set(exclude_target_ids or set())
    sites = []
    for seq_id, raw_seq in reference_sequences.items():
        if seq_id in excluded:
            continue
        seq = to_dna(raw_seq)
        for possible in _scan_spcas9_sites(seq, pam=pam):
            item = dict(possible)
            item["target_id"] = seq_id
            item["reference_length"] = len(seq)
            sites.append(item)
    return SgRNAReferenceSites(sites)


def _candidate_reference_locus(candidate: Dict) -> Dict[str, int]:
    position_start = int(candidate.get("source_position_start", candidate.get("position_start", 0)) or 0)
    position_end = int(candidate.get("source_position_end", candidate.get("position_end", position_start)) or position_start)
    pam_start = int(candidate.get("source_pam_start", candidate.get("pam_start", position_end)) or position_end)
    pam_end = int(candidate.get("source_pam_end", candidate.get("pam_end", pam_start)) or pam_start)
    return {
        "position_start": position_start,
        "position_end": position_end,
        "pam_start": pam_start,
        "pam_end": pam_end,
        "locus_start": min(position_start, pam_start),
        "locus_end": max(position_end, pam_end),
    }


def score_sgrna_offtargets_from_sites(
    candidate: Dict,
    reference_sites,
    exclude_target_id: Optional[str] = None,
    exclude_target_ids: Optional[Iterable[str]] = None,
    max_mismatches: int = 5,
) -> Dict:
    """Score a guide against precomputed SpCas9 NRG reference sites."""
    spacer = to_dna(candidate["spacer_dna"])
    seed_12 = spacer[-12:]
    intended = _candidate_reference_locus(candidate)
    intended_strand = candidate.get("strand")
    hits = []
    fully_excluded_ids = set(exclude_target_ids or set())
    if exclude_target_id in fully_excluded_ids:
        fully_excluded_ids.remove(exclude_target_id)
    if hasattr(reference_sites, "candidate_sites"):
        possible_sites = reference_sites.candidate_sites(spacer, max_mismatches)
    else:
        possible_sites = reference_sites
    for possible in possible_sites:
        seq_id = possible["target_id"]
        if seq_id in fully_excluded_ids:
            continue
        possible_locus_start = min(possible["position_start"], possible["pam_start"])
        possible_locus_end = max(possible["position_end"], possible["pam_end"])
        overlaps_intended_locus = (
            seq_id == exclude_target_id
            and max(intended["locus_start"], possible_locus_start) < min(intended["locus_end"], possible_locus_end)
        )
        is_intended_site = (
            seq_id == exclude_target_id
            and possible.get("position_start") == intended["position_start"]
            and possible.get("strand") == intended_strand
        )
        if is_intended_site or overlaps_intended_locus:
            continue
        target_spacer = possible["spacer_dna"]
        mismatches = _mismatch_positions(spacer, target_spacer)
        if len(mismatches) > max_mismatches:
            continue
        cfd = _cfd_like_score(mismatches)
        seed_12_match = target_spacer[-12:] == seed_12
        mismatch_floor = {0: 100.0, 1: 72.0, 2: 45.0, 3: 24.0, 4: 8.0, 5: 4.0}.get(len(mismatches), 0.0)
        risk_score = max(cfd, mismatch_floor)
        if seed_12_match:
            risk_score += 18
        if possible["pam"][1] == "A":
            risk_score *= 0.85
        risk_score = round(risk_score, 2)
        reasons = [f"{len(mismatches)} mismatch Cas9 site", f"PAM {possible['pam']}"]
        if possible["pam"][1] == "A":
            reasons.append("alternative NRG/NAG PAM")
        if seed_12_match:
            reasons.append("12 nt PAM-proximal seed match")
        if any(pos >= 13 for pos in mismatches):
            reasons.append("PAM-proximal mismatch")
        if not mismatches:
            reasons.append("perfect spacer match")
        reference_length = int(possible.get("reference_length", 0) or 0)
        validation_start = max(0, min(possible["position_start"], possible["pam_start"]) - 100)
        validation_end = min(reference_length, max(possible["position_end"], possible["pam_end"]) + 100)
        hits.append({
            "target_id": seq_id,
            "risk_score": min(100.0, risk_score),
            "cfd_like_score": cfd,
            "mismatches": len(mismatches),
            "mismatch_positions": mismatches,
            "pam": possible["pam"],
            "genomic_pam": possible.get("genomic_pam", possible["pam"]),
            "strand": possible["strand"],
            "position": possible["position_start"],
            "cut_site": possible["cut_site"],
            "target_spacer": target_spacer,
            "target_protospacer_pam": target_spacer + possible["pam"],
            "seed_12_match": seed_12_match,
            "validation_window_start": validation_start,
            "validation_window_end": validation_end,
            "reasons": reasons,
        })
    return _risk_from_hits(hits, max_mismatches=max_mismatches)


def score_sgrna_offtargets(
    candidate: Dict,
    reference_sequences: Dict[str, str],
    exclude_target_id: Optional[str] = None,
    exclude_target_ids: Optional[Iterable[str]] = None,
    max_mismatches: int = 5,
) -> Dict:
    """Search SpCas9 NRG-adjacent off-targets with <= max_mismatches mismatches."""
    fully_excluded_ids = set(exclude_target_ids or set())
    if exclude_target_id in fully_excluded_ids:
        fully_excluded_ids.remove(exclude_target_id)
    reference_sites = build_sgrna_reference_sites(reference_sequences, exclude_target_ids=fully_excluded_ids)
    return score_sgrna_offtargets_from_sites(
        candidate,
        reference_sites,
        exclude_target_id=exclude_target_id,
        max_mismatches=max_mismatches,
    )


def _restriction_site_warnings(spacer: str) -> List[Dict]:
    warnings = []
    for enzyme, sites in RESTRICTION_ENZYME_SITES.items():
        for site in sites:
            pos = spacer.find(site)
            if pos >= 0:
                warnings.append({
                    "enzyme": enzyme,
                    "site": site,
                    "position": pos + 1,
                    "message": f"{enzyme} restriction enzyme site {site} occurs in spacer at position {pos + 1}",
                })
    return warnings


def design_sgrna_cloning_oligos(
    spacer_dna: str,
    forward_overhang: str = "CACC",
    reverse_overhang: str = "AAAC",
) -> Dict:
    """Return common BbsI-compatible sgRNA cloning oligos."""
    spacer = to_dna(spacer_dna)
    if not _valid_spacer(spacer):
        raise ValueError("sgRNA spacer must be exactly 20 A/C/G/T bases")
    u6_spacer = spacer if spacer.startswith("G") else "G" + spacer
    reverse = reverse_complement_dna(u6_spacer)
    custom_reverse = reverse_complement_dna(spacer)
    warnings = _restriction_site_warnings(spacer)
    notes = [
        "BbsI-style oligos compatible with common SpCas9 sgRNA cloning vectors such as pX330/lentiCRISPR-style workflows; verify vector overhangs before ordering."
    ]
    if warnings:
        notes.append(
            "WARNING: spacer contains a restriction enzyme site used by common sgRNA cloning workflows; verify the selected vector and cloning enzyme before ordering."
        )
    return {
        "spacer_dna": spacer,
        "guide_rna": spacer.replace("T", "U"),
        "u6_spacer_dna": u6_spacer,
        "px330_forward_oligo": {
            "sequence": "CACC" + u6_spacer,
            "role": "BbsI forward cloning oligo",
            "length": len("CACC" + u6_spacer),
        },
        "px330_reverse_oligo": {
            "sequence": "AAAC" + reverse + "C",
            "role": "BbsI reverse cloning oligo",
            "length": len("AAAC" + reverse + "C"),
        },
        "custom_forward_oligo": {
            "sequence": to_dna(forward_overhang) + spacer,
            "role": "Custom forward cloning oligo",
            "length": len(to_dna(forward_overhang) + spacer),
        },
        "custom_reverse_oligo": {
            "sequence": to_dna(reverse_overhang) + custom_reverse,
            "role": "Custom reverse cloning oligo",
            "length": len(to_dna(reverse_overhang) + custom_reverse),
        },
        "restriction_site_warnings": warnings,
        "notes": " ".join(notes),
    }


def design_genotyping_primers(sequence: str, cut_site: int) -> Dict:
    """Suggest simple PCR primers around the Cas9 cut site when enough context exists."""
    dna = to_dna(sequence)
    cut_site = max(0, min(len(dna), int(cut_site)))
    start = max(0, cut_site - 120)
    end = min(len(dna), cut_site + 120)
    amplicon = dna[start:end]
    if len(amplicon) < 80:
        return {}
    try:
        primer_set = design_t7_primers(amplicon, result_id="sgRNA_genotyping")
    except ValueError:
        return {}
    return {
        "amplicon_start": start,
        "amplicon_end": end,
        "amplicon_length": len(amplicon),
        "forward_primer": primer_set["forward_primer"],
        "reverse_primer": primer_set["reverse_primer"],
        "notes": "PCR primers flank the predicted Cas9 cut site for genotyping; validate specificity against the full genome before ordering.",
    }
