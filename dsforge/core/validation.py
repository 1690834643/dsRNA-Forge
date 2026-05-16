"""
Build experiment-oriented off-target validation records.
"""

from typing import Dict, List

from dsforge.core.sequence import normalize_sequence


def _longest_contiguous_match(query: str, subject: str) -> int:
    """Return the longest exact contiguous match between query and subject."""
    query = normalize_sequence(query)
    subject = normalize_sequence(subject)
    if not query or not subject:
        return 0
    previous = [0] * (len(subject) + 1)
    best = 0
    for q_base in query:
        current = [0] * (len(subject) + 1)
        for j, s_base in enumerate(subject, start=1):
            if q_base == s_base:
                current[j] = previous[j - 1] + 1
                best = max(best, current[j])
        previous = current
    return best


def _mismatch_count(query: str, subject: str) -> int:
    length = min(len(query), len(subject))
    if length == 0:
        return max(len(query), len(subject))
    return sum(a != b for a, b in zip(query[:length], subject[:length])) + abs(len(query) - len(subject))


def _best_match_for_target(matches: List[Dict], target_id: str, position=None) -> Dict:
    candidates = [m for m in matches if m.get("target_id") == target_id]
    if position is not None:
        positioned = [m for m in candidates if int(m.get("position", -1) or -1) == int(position)]
        if positioned:
            candidates = positioned
    if not candidates:
        return {"target_id": target_id, "match_type": "ranked_risk", "length": 0, "position": 0}
    return max(
        candidates,
        key=lambda item: (
            str(item.get("match_type", "")).startswith("Cas9_seed12_POT"),
            int(item.get("length", 0) or 0),
        ),
    )


def build_validation_hits(result: Dict, transcriptome, max_hits: int = 5) -> List[Dict]:
    """Return ranked off-target snippets and validation directions for one result."""
    off_target = result.get("off_target") or {}
    top_targets = off_target.get("top_targets") or []
    matches = off_target.get("matches") or []
    query = normalize_sequence(result.get("sequence", result.get("candidate_seq", "")))
    if not query:
        return []

    hits = []
    for target in top_targets[:max_hits]:
        target_id = target.get("target_id", "")
        subject = transcriptome.get_sequence(target_id) if transcriptome is not None else ""
        if not subject:
            continue
        subject = normalize_sequence(subject)
        match = _best_match_for_target(matches, target_id, position=target.get("position"))
        position = int(match.get("position", 0) or 0)
        window_start_value = match.get("validation_window_start")
        window_end_value = match.get("validation_window_end")
        window_start = int(position if window_start_value is None else window_start_value)
        window_end = int((window_start + len(query)) if window_end_value is None else window_end_value)
        window_start = max(0, min(window_start, len(subject)))
        window_end = max(window_start, min(window_end, len(subject)))
        target_fragment = subject[window_start:window_end]
        if len(target_fragment) < min(len(query), 8):
            window_start = max(0, position - 10)
            window_end = min(len(subject), position + len(query) + 10)
            target_fragment = subject[window_start:window_end]

        seed = query[1:8] if len(query) >= 8 else query
        matched_spacer = normalize_sequence(match.get("target_spacer", "")) or subject[position:position + len(query)]
        longest = _longest_contiguous_match(query, target_fragment)
        mismatch = _mismatch_count(query, matched_spacer[: len(query)])
        match_type = match.get("match_type", "ranked_risk")
        action = "qPCR + local alignment"
        if match_type.startswith("Cas9_seed12_POT"):
            action = "Priority PCR amplicon sequencing for seed-matched Cas9 POT site; inspect with ICE/TIDE/Sanger"
        elif match_type.startswith("Cas9_"):
            action = "Targeted PCR amplicon sequencing; inspect editing with ICE/TIDE/Sanger"
        elif "20bp" in match_type or longest >= 20:
            action = "High priority: avoid or validate by BLAST/qPCR"
        elif "16bp" in match_type or longest >= 16:
            action = "Validate by local alignment and qPCR"
        elif "seed" in match_type or seed in target_fragment:
            action = "Check seed-mediated off-target risk"

        hits.append({
            "target_id": target_id,
            "risk_score": target.get("risk_score", 0),
            "reasons": target.get("reasons", []),
            "match_type": match_type,
            "target_position": position,
            "query_fragment": query[: len(target_fragment)] if len(target_fragment) < len(query) else query,
            "target_fragment": target_fragment,
            "matched_spacer": matched_spacer,
            "target_protospacer_pam": match.get("target_protospacer_pam", ""),
            "validation_window_start": window_start,
            "validation_window_end": window_end,
            "seed": seed,
            "seed_match": seed in target_fragment,
            "longest_contiguous_match": longest,
            "mismatch_count": mismatch,
            "validation_action": action,
        })

    hits.sort(key=lambda item: (item["risk_score"], item["longest_contiguous_match"]), reverse=True)
    return hits
