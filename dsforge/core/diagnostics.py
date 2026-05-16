"""
User-facing diagnosis for empty or over-filtered design results.
"""

from typing import Dict, List

from dsforge.core.sequence import AMBIGUOUS_RNA_BASES, generate_candidates, normalize_sequence


def _window_for_mode(mode: str, config) -> tuple[int, int]:
    if mode == "siRNA":
        return 21, 21
    if mode == "DsiRNA":
        return 27, 27
    if mode == "sgRNA":
        return 23, 23
    return getattr(config, "length_min", 21), getattr(config, "length_max", 21)


def _diagnose_empty_sgrna(seq: str) -> Dict[str, List[str]]:
    dna = seq.replace("U", "T")
    valid_bases = set("ACGT")
    has_plus_pam = any(dna[i + 1 : i + 3] == "GG" for i in range(0, max(0, len(dna) - 2)))
    has_reverse_pam = any(dna[i : i + 2] == "CC" for i in range(0, max(0, len(dna) - 2)))

    if not has_plus_pam and not has_reverse_pam:
        return {
            "reasons": ["目标序列中没有可用 SpCas9 PAM：正向需要 NGG，反向在输入正链上表现为 CCN。"],
            "suggestions": ["换一个包含 NGG/CCN 的目标区域，或上传带上下游 flanking sequence 的更长片段。"],
        }

    if set(dna) - valid_bases:
        return {
            "reasons": ["检测到 SpCas9 PAM，但 PAM 附近的 20 nt spacer 含有 N/模糊碱基，无法生成精确 sgRNA。"],
            "suggestions": ["使用 A/C/G/T 明确碱基的参考序列，或避开含 N 的区域。"],
        }

    return {
        "reasons": ["检测到 SpCas9 PAM，但没有足够的 20 nt spacer 上下文可生成候选。"],
        "suggestions": ["上传更长的目标片段，确保 NGG 前方或 CCN 后方至少有 20 nt 可用序列。"],
    }


def diagnose_design_outcome(
    target_seq: str,
    config,
    total_candidates: int,
    passed_candidates: int,
) -> Dict[str, List[str]]:
    """Return plain Chinese reasons and suggestions for weak design outcomes."""
    seq = normalize_sequence(target_seq)
    reasons: List[str] = []
    suggestions: List[str] = []

    mode = getattr(config, "mode", "siRNA")
    min_len, max_len = _window_for_mode(mode, config)

    if len(seq) < min_len:
        if mode == "sgRNA":
            reasons.append(f"目标序列太短：当前长度 {len(seq)} nt，sgRNA 模式至少需要 23 nt（20 nt spacer + 3 nt PAM）。")
        else:
            reasons.append(f"目标序列太短：当前长度 {len(seq)} nt，{mode} 模式至少需要 {min_len} nt。")
        if mode == "long_dsRNA":
            suggestions.append("换用 siRNA/DsiRNA 模式，或提供更长的目标序列。")
        elif mode == "sgRNA":
            suggestions.append("提供包含目标位点上下游序列的片段，并确保其中有 NGG/CCN PAM。")
        else:
            suggestions.append("提供更长的目标序列，或检查粘贴/上传的 FASTA 是否完整。")
        return {"reasons": reasons, "suggestions": suggestions}

    if total_candidates == 0:
        if mode == "sgRNA":
            return _diagnose_empty_sgrna(seq)
        ambiguous = sorted(set(seq) & AMBIGUOUS_RNA_BASES)
        unconstrained = list(generate_candidates(seq, mode, min_len, max_len, 0, 100, exclude_poly=9999, max_candidates=100000))
        gc_only = list(
            generate_candidates(
                seq,
                mode,
                min_len,
                max_len,
                getattr(config, "gc_min", 0),
                getattr(config, "gc_max", 100),
                exclude_poly=9999,
                max_candidates=100000,
            )
        )
        if not unconstrained and ambiguous:
            reasons.append(
                f"目标序列包含 N/模糊 IUPAC 碱基（{', '.join(ambiguous)}）；含这些碱基的候选窗口会被保守跳过。"
            )
            suggestions.append("换用 A/U/G/C 明确的参考序列，或上传避开 N/模糊碱基的目标片段。")
        elif not unconstrained:
            reasons.append("目标序列长度不满足当前设计模式的候选窗口要求。")
            suggestions.append("切换设计模式，或在高级设置里调整长度范围。")
        elif not gc_only:
            reasons.append("GC 含量范围过窄，所有候选都被过滤。")
            suggestions.append("在高级设置里放宽 GC 范围，例如 20-80%。")
        else:
            reasons.append("连续重复碱基过滤过严，所有候选都被过滤。")
            suggestions.append("在高级设置里放宽过滤条件，或检查目标序列是否低复杂度。")
    elif passed_candidates == 0:
        reasons.append("已生成候选，但没有候选通过当前评分或脱靶筛查。")
        suggestions.append("优先查看候选表中的分数；可尝试放宽 GC 范围或换一个目标区域。")
        suggestions.append("如果转录组里存在高度相似基因，脱靶筛查可能会过滤掉全部候选。")

    if not reasons:
        reasons.append("设计完成，但当前结果较少。")
        suggestions.append("可导出结果后优先选择排名靠前且通过筛选的候选。")

    return {"reasons": reasons, "suggestions": suggestions}
