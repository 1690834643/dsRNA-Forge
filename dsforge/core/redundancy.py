"""
Collapse highly overlapping sliding-window candidates into representative recommendations.
"""

from typing import Dict, List


def interval_jaccard(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    """Return intersection-over-union for two half-open intervals."""
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    if overlap <= 0:
        return 0.0
    union = max(a_end, b_end) - min(a_start, b_start)
    return overlap / union if union else 0.0


def _score(result: Dict) -> float:
    return float(result.get("recommendation_score", result.get("consensus_score", 0)) or 0)


def cluster_redundant_results(results: List[Dict], overlap_threshold: float = 0.8) -> List[Dict]:
    """
    Greedily keep the best-scoring representative for each highly overlapping interval cluster.

    Input order is not trusted; results are sorted by recommendation score first. A shifted
    1 bp siRNA/long-dsRNA window clusters with its neighbor, while a much longer nested
    interval remains separate because Jaccard overlap is low.
    """
    ranked = sorted(results, key=_score, reverse=True)
    representatives: List[Dict] = []
    cluster_members: List[List[Dict]] = []

    for result in ranked:
        start = int(result.get("position_start", 0) or 0)
        end = int(result.get("position_end", start) or start)
        assigned = False
        for idx, rep in enumerate(representatives):
            rep_start = int(rep.get("position_start", 0) or 0)
            rep_end = int(rep.get("position_end", rep_start) or rep_start)
            if interval_jaccard(start, end, rep_start, rep_end) >= overlap_threshold:
                cluster_members[idx].append(result)
                assigned = True
                break
        if not assigned:
            representatives.append(dict(result))
            cluster_members.append([result])

    for cluster_id, (rep, members) in enumerate(zip(representatives, cluster_members), start=1):
        starts = [int(m.get("position_start", 0) or 0) for m in members]
        ends = [int(m.get("position_end", 0) or 0) for m in members]
        rep["cluster_id"] = cluster_id
        rep["cluster_size"] = len(members)
        rep["alternative_count"] = max(0, len(members) - 1)
        rep["cluster_span"] = f"{min(starts)}-{max(ends)}" if starts and ends else rep.get("position", "")

    representatives.sort(key=_score, reverse=True)
    for rank, rep in enumerate(representatives, start=1):
        rep["rank"] = rank
    return representatives
