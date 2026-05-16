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
RESTRICTION_ENZYME_SITES = {
    "BsaI": ("GGTCTC", "GAGACC"),
    "BbsI": ("GAAGAC", "GTCTTC"),
    "BsmBI": ("CGTCTC", "GAGACG"),
}


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


def score_sgrna_offtargets(
    candidate: Dict,
    reference_sequences: Dict[str, str],
    exclude_target_id: Optional[str] = None,
    max_mismatches: int = 5,
) -> Dict:
    """Search SpCas9 NRG-adjacent off-targets with <= max_mismatches mismatches."""
    spacer = to_dna(candidate["spacer_dna"])
    seed_12 = spacer[-12:]
    intended_start = candidate.get("position_start")
    intended_strand = candidate.get("strand")
    intended_locus_start = min(
        int(candidate.get("position_start", 0) or 0),
        int(candidate.get("pam_start", candidate.get("position_start", 0)) or 0),
    )
    intended_locus_end = max(
        int(candidate.get("position_end", intended_locus_start) or intended_locus_start),
        int(candidate.get("pam_end", candidate.get("position_end", intended_locus_start)) or intended_locus_start),
    )
    hits = []
    for seq_id, raw_seq in reference_sequences.items():
        seq = to_dna(raw_seq)
        for possible in _scan_spcas9_sites(seq, pam="NRG"):
            possible_locus_start = min(possible["position_start"], possible["pam_start"])
            possible_locus_end = max(possible["position_end"], possible["pam_end"])
            overlaps_intended_locus = (
                seq_id == exclude_target_id
                and max(intended_locus_start, possible_locus_start) < min(intended_locus_end, possible_locus_end)
            )
            is_intended_site = (
                seq_id == exclude_target_id
                and possible.get("position_start") == intended_start
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
            validation_start = max(0, min(possible["position_start"], possible["pam_start"]) - 100)
            validation_end = min(len(seq), max(possible["position_end"], possible["pam_end"]) + 100)
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
