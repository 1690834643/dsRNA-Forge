"""
Simple T7 primer design for long dsRNA wet-lab handoff.
"""

import math
from typing import Dict, List, Optional, Set, Tuple

from dsforge.core.sequence import gc_content, normalize_sequence


def check_primer_specificity(
    primer_dna: str,
    transcriptome_sequences: Dict[str, str],
    max_mismatches: int = 2,
    max_matches_to_report: int = 5,
    exclude_ids: Optional[set] = None,
    three_prime_anchor: int = 8,
    three_prime_risk_len: int = 10,
) -> Dict:
    """Check primer off-target binding on both transcript strands.

    Args:
        primer_dna: Primer sequence (DNA)
        transcriptome_sequences: Dict of seq_id -> sequence
        max_mismatches: Maximum tolerated mismatches for a full-length off-target call
        max_matches_to_report: Cap on reported matches for performance
        exclude_ids: Set of seq_ids to exclude (e.g. the target gene itself)
        three_prime_anchor: Exact 3' suffix length used to pre-index likely annealing sites
        three_prime_risk_len: 3' exact-match length that can rescue otherwise borderline sites
    """
    primer_dna = _to_dna(primer_dna)
    if not primer_dna:
        return {"specific": True, "matches": [], "warning": ""}

    excluded = exclude_ids or set()
    min_anchor = max(4, min(three_prime_anchor, len(primer_dna)))
    max_borderline_mismatches = max_mismatches + 2
    orientations: Tuple[Tuple[str, str], ...] = (
        ("+", primer_dna),
        ("-", _reverse_complement_dna(primer_dna)),
    )
    matches: List[Dict] = []
    seen: Set[Tuple[str, int, str]] = set()

    for seq_id, seq in transcriptome_sequences.items():
        if seq_id in excluded:
            continue
        seq = seq.upper().replace("U", "T")
        if len(seq) < len(primer_dna):
            continue
        for strand, query in orientations:
            for anchor, anchor_offset in _primer_anchor_specs(query, max_mismatches, min_anchor):
                start = seq.find(anchor)
                while start != -1:
                    i = start - anchor_offset
                    if i < 0 or i + len(query) > len(seq):
                        start = seq.find(anchor, start + 1)
                        continue

                    key = (seq_id, i, strand)
                    if key in seen:
                        start = seq.find(anchor, start + 1)
                        continue
                    seen.add(key)

                    window = seq[i : i + len(query)]
                    if set(window) - set("ACGT"):
                        start = seq.find(anchor, start + 1)
                        continue
                    mismatches = sum(a != b for a, b in zip(query, window))
                    three_prime_match_len = _three_prime_match_len(query, window)
                    is_high_risk = (
                        mismatches <= max_mismatches
                        or (
                            three_prime_match_len >= three_prime_risk_len
                            and mismatches <= max_borderline_mismatches
                        )
                    )
                    if not is_high_risk:
                        start = seq.find(anchor, start + 1)
                        continue
                    risk_score = _primer_offtarget_risk_score(mismatches, three_prime_match_len)
                    matches.append({
                        "target_id": seq_id,
                        "position": i,
                        "mismatches": mismatches,
                        "window": window,
                        "strand": strand,
                        "three_prime_match_len": three_prime_match_len,
                        "risk_score": risk_score,
                        "reason": _primer_offtarget_reason(mismatches, three_prime_match_len, max_mismatches),
                    })
                    start = seq.find(anchor, start + 1)

    matches.sort(
        key=lambda item: (
            -item["risk_score"],
            item["mismatches"],
            -item["three_prime_match_len"],
            item["target_id"],
            item["position"],
            item["strand"],
        )
    )
    reported = matches[:max_matches_to_report]
    return {
        "specific": len(matches) == 0,
        "matches": reported,
        "total_matches": len(matches),
        "parameters": {
            "max_mismatches": max_mismatches,
            "three_prime_anchor": three_prime_anchor,
            "three_prime_risk_len": three_prime_risk_len,
        },
        "warning": (
            "Non-specific: potential off-target amplification sites detected. "
            "Verify with local BLAST or NCBI Primer-BLAST before ordering."
        ) if matches else "",
    }


T7_PROMOTER = "TAATACGACTCACTATAGGG"
DNA_COMPLEMENT = str.maketrans("AUGCT", "TACGA")


def _to_dna(sequence: str) -> str:
    return normalize_sequence(sequence).replace("U", "T")


def _reverse_complement_dna(sequence: str) -> str:
    return _to_dna(sequence).translate(DNA_COMPLEMENT)[::-1]


def _tm_wallace(sequence: str) -> float:
    sequence = _to_dna(sequence)
    return float(2 * (sequence.count("A") + sequence.count("T")) + 4 * (sequence.count("G") + sequence.count("C")))


_NN_THERMODYNAMICS = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}


def calculate_primer_tm(
    sequence: str,
    sodium_molar: float = 0.05,
    primer_molar: float = 250e-9,
) -> float:
    """Return an approximate DNA nearest-neighbor melting temperature in Celsius."""
    dna = _to_dna(sequence)
    if len(dna) < 8 or set(dna) - set("ACGT"):
        return round(_tm_wallace(dna), 1)

    dh = 0.2
    ds = -5.7
    for i in range(len(dna) - 1):
        pair = dna[i : i + 2]
        if pair not in _NN_THERMODYNAMICS:
            return round(_tm_wallace(dna), 1)
        pair_dh, pair_ds = _NN_THERMODYNAMICS[pair]
        dh += pair_dh
        ds += pair_ds

    salt = max(float(sodium_molar), 1e-6)
    concentration = max(float(primer_molar), 1e-12)
    tm = (dh * 1000.0) / (ds + 1.987 * math.log(concentration / 4.0))
    tm = tm - 273.15 + 16.6 * math.log10(salt)
    return round(tm, 1)


def _three_prime_match_len(query: str, window: str) -> int:
    count = 0
    for q_base, w_base in zip(reversed(query), reversed(window)):
        if q_base != w_base:
            break
        count += 1
    return count


def _primer_anchor_specs(query: str, max_mismatches: int, three_prime_anchor: int) -> List[Tuple[str, int]]:
    """Return exact anchors that guarantee discovery of <=max_mismatches near matches."""
    specs: List[Tuple[str, int]] = []
    seen = set()

    def add(start: int, end: int):
        anchor = query[start:end]
        key = (anchor, start)
        if len(anchor) < 4 or key in seen:
            return
        seen.add(key)
        specs.append((anchor, start))

    add(len(query) - three_prime_anchor, len(query))
    segments = max(1, int(max_mismatches) + 1)
    for segment in range(segments):
        start = round(segment * len(query) / segments)
        end = round((segment + 1) * len(query) / segments)
        add(start, end)
    return specs


def _primer_offtarget_risk_score(mismatches: int, three_prime_match_len: int) -> int:
    score = 100 - mismatches * 18 + min(30, three_prime_match_len * 2)
    return max(0, min(100, score))


def _primer_offtarget_reason(mismatches: int, three_prime_match_len: int, max_mismatches: int) -> str:
    if mismatches <= max_mismatches:
        return "full-length near match"
    return "strong 3' anchored match"


def _primer_record(sequence: str, role: str) -> Dict:
    sequence = _to_dna(sequence)
    return {
        "role": role,
        "sequence": sequence,
        "length": len(sequence),
        "tm": calculate_primer_tm(sequence),
        "gc_percent": round(gc_content(sequence), 1),
    }


def _pick_forward(sequence: str) -> str:
    dna = _to_dna(sequence)
    for length in range(20, 25):
        primer = dna[:length]
        gc = gc_content(primer)
        tm = calculate_primer_tm(primer)
        if 35 <= gc <= 65 and 50 <= tm <= 72:
            return primer
    return dna[:22] if len(dna) >= 22 else dna


def _pick_reverse(sequence: str) -> str:
    dna = _to_dna(sequence)
    for length in range(20, 25):
        primer = _reverse_complement_dna(dna[-length:])
        gc = gc_content(primer)
        tm = calculate_primer_tm(primer)
        if 35 <= gc <= 65 and 50 <= tm <= 72:
            return primer
    return _reverse_complement_dna(dna[-22:] if len(dna) >= 22 else dna)


def design_t7_primers(sequence: str, result_id: Optional[str] = None) -> Dict:
    """Design basic PCR and T7-promoter primers for one dsRNA candidate."""
    dna = _to_dna(sequence)
    if len(dna) < 40:
        raise ValueError("T7 primer design requires at least 40 nt of candidate sequence")

    forward = _pick_forward(dna)
    reverse = _pick_reverse(dna)
    return {
        "result_id": result_id or "",
        "product_length": len(dna),
        "forward_primer": _primer_record(forward, "forward"),
        "reverse_primer": _primer_record(reverse, "reverse"),
        "t7_forward_primer": _primer_record(T7_PROMOTER + forward, "T7 forward"),
        "t7_reverse_primer": _primer_record(T7_PROMOTER + reverse, "T7 reverse"),
        "notes": "Use T7-forward + T7-reverse for in vitro transcription templates; verify amplicon uniqueness before ordering.",
    }
