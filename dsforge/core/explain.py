"""
Result explanation helpers for experiment-facing recommendations.
"""

from typing import Dict, Iterable, List, Optional


def _risk_label(off_target: Dict) -> str:
    level = (off_target or {}).get("risk_level", "low")
    score = (off_target or {}).get("risk_score", 0)
    return f"{str(level).upper()} risk, score {score}"


def explain_result(result: Dict, mode: str = "") -> Dict:
    """Create a structured plain-language explanation for one recommendation."""
    off_target = result.get("off_target") or {}
    score = float(result.get("consensus_score", 0) or 0)
    recommendation = float(result.get("recommendation_score", score) or score)
    cluster_size = int(result.get("cluster_size", 1) or 1)
    passed = bool(result.get("passed", False))
    top_targets = off_target.get("top_targets") or []

    efficacy_notes = [
        f"Consensus score {score:.1f}; ranked recommendation score {recommendation:.1f}.",
    ]
    if cluster_size > 1:
        efficacy_notes.append(f"Cluster Size {cluster_size}: represents adjacent similar windows.")
    if mode == "long_dsRNA" and result.get("pool"):
        pool = result["pool"]
        efficacy_notes.append(
            f"Dicer pool score {pool.get('pool_score', score):.1f}; "
            f"{len(pool.get('product_details', []))} heuristic siRNA products for relative ranking."
        )
    if mode == "sgRNA" and result.get("sgrna"):
        sg = result["sgrna"]
        genomic_pam = sg.get("genomic_pam", sg.get("pam", ""))
        genomic_text = f", genomic PAM {genomic_pam}" if genomic_pam and genomic_pam != sg.get("pam", "") else ""
        efficacy_notes.append(
            f"SpCas9 guide: spacer {sg.get('spacer_dna', '')}, PAM {sg.get('pam', '')}, "
            f"strand {sg.get('strand', '')}{genomic_text}, cut site {sg.get('cut_site', '')}."
        )

    risk_notes: List[str] = [f"Off-target summary: {_risk_label(off_target)}."]
    summary = off_target.get("summary") or {}
    mismatch_counts = summary.get("mismatch_counts") or {}
    pot_counts = summary.get("pot_mismatch_counts") or {}
    if mismatch_counts:
        risk_notes.append(
            "sgRNA OT 0M-5M counts: "
            + ", ".join(f"{bucket}={mismatch_counts.get(bucket, 0)}" for bucket in ["0M", "1M", "2M", "3M", "4M", "5M"])
            + f"; NRG/NAG PAM hits={summary.get('nrg_pam_hits', 0)}."
        )
    if pot_counts:
        risk_notes.append(
            "sgRNA seed12 POT 0M-5M counts: "
            + ", ".join(f"{bucket}={pot_counts.get(bucket, 0)}" for bucket in ["0M", "1M", "2M", "3M", "4M", "5M"])
            + f"; total POT={summary.get('total_pot', 0)}."
        )
    for target in top_targets[:5]:
        reasons = ", ".join(target.get("reasons", []))
        suffix = f" ({reasons})" if reasons else ""
        risk_notes.append(
            f"Top risk target {target.get('target_id', '')}: "
            f"score {target.get('risk_score', '')}{suffix}."
        )
    if off_target.get("validation_direction"):
        risk_notes.append(off_target["validation_direction"])

    method_notes = []
    rnaup = result.get("rnaup") or {}
    rnaup_method = (rnaup.get("details") or {}).get("method", "")
    if rnaup_method == "RNAup-cli":
        method_notes.append(f"RNAup CLI refinement available; ΔG={rnaup.get('dg')}.")
    elif rnaup_method:
        method_notes.append(
            f"RNAduplex fallback used for RNAup refinement; ΔG={rnaup.get('dg')}. "
            "Treat this as screening evidence, not full RNAup precision."
        )
    elif result.get("thermodynamics"):
        method_notes.append("RNAduplex/RNAcofold thermodynamic screening was used.")
    else:
        method_notes.append("Thermodynamic refinement was not available for this candidate.")
    if mode == "sgRNA":
        method_notes.append(
            "sgRNA off-target scan is limited to currently loaded reference/background sequences; "
            "load genome FASTA as an extra background for genome-scale Cas9 off-target screening."
        )

    if passed:
        decision = "优先候选" if off_target.get("risk_level", "low") == "low" else "可验证候选"
    else:
        decision = "谨慎使用或替换候选"

    validation_notes = [
        "优先对目标基因 knockdown 效率做 qPCR 验证。",
        "对 Top risk targets 做 BLAST/局部比对复核；中高风险对象建议追加 qPCR。",
    ]
    if mode == "long_dsRNA":
        validation_notes.append("长 dsRNA 的 Dicer pool 为简化模拟；建议先确认扩增片段唯一性、T7 体外转录产物长度和 knockdown 效率。")
    if mode == "sgRNA":
        validation_notes.append("sgRNA 建议对 Top off-target 位点做 PCR 扩增测序，并用 ICE/TIDE 或 Sanger 分析确认编辑。")
        validation_notes.append("如果当前只加载了转录组，Cas9 脱靶报告不覆盖 intron/intergenic 区域；需要基因组 FASTA 复核。")

    return {
        "decision": decision,
        "summary": f"推荐理由：{decision}；{_risk_label(off_target)}；score {recommendation:.1f}.",
        "efficacy_notes": efficacy_notes,
        "risk_notes": risk_notes,
        "method_notes": method_notes,
        "validation_notes": validation_notes,
    }


def render_region_map(
    target_length: int,
    start: int,
    end: int,
    risk_positions: Optional[Iterable[int]] = None,
    width: int = 60,
) -> str:
    """Render a compact ASCII target-region map for the result detail panel."""
    target_length = max(1, int(target_length or 1))
    start = max(0, min(int(start or 0), target_length))
    end = max(start, min(int(end or start), target_length))
    width = max(20, int(width or 60))
    risk_positions = list(risk_positions or [])

    chars = ["-"] * width
    map_start = int(start / target_length * width)
    map_end = max(map_start + 1, int(end / target_length * width))
    for i in range(map_start, min(width, map_end)):
        chars[i] = "#"
    for pos in risk_positions:
        idx = max(0, min(width - 1, int(pos / target_length * width)))
        chars[idx] = "!"

    return (
        f"target 0 {'-' * max(1, width - 12)} {target_length} nt\n"
        f"candidate [{''.join(chars)}]\n"
        "# candidate region, ! off-target/risk signal"
    )
