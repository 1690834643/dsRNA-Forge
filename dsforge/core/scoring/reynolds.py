"""
Reynolds 评分器 (2004)
基于 180 条 siRNA 的 8 条评分标准
"""

from dsforge.core.scoring.base import ScoringRule, register_rule


@register_rule
class ReynoldsRule(ScoringRule):
    """
    Reynolds et al. (2004) — Rational siRNA design for RNA interference
    """

    name = "reynolds"
    description = "Reynolds rules (2004): 8 criteria, cutoff >= 6"
    CUTOFF = 6

    def score(self, sequence: str) -> dict:
        """
        对 21nt siRNA antisense 序列评分
        位置 1 = antisense 5' 端, 位置 19 = antisense 3' 端 (靠近切割位点)
        """
        if not self.validate_sequence(sequence):
            return {"score": 0.0, "passed": False, "violations": ["Invalid sequence"], "details": {}}

        seq = sequence.upper().replace("T", "U")
        if len(seq) < 19:
            return {"score": 0.0, "passed": False, "violations": ["Sequence too short (< 19nt)"], "details": {}}

        # 只取前 19nt 进行评分（经典规则针对 19nt 核心区域）
        core = seq[:19]
        score = 0
        violations = []
        details = {}

        # 1. GC 含量 30%–52%
        gc_count = core.count("G") + core.count("C")
        gc_pct = (gc_count / 19) * 100
        if 30 <= gc_pct <= 52:
            score += 1
            details["gc_content"] = {"value": gc_pct, "passed": True}
        else:
            violations.append(f"GC content {gc_pct:.1f}% not in 30-52%")
            details["gc_content"] = {"value": gc_pct, "passed": False}

        # 2. 15–19 位至少 3 个 A/U
        tail = core[14:19]  # 位置 15-19 (0-indexed 14-18)
        au_count = tail.count("A") + tail.count("U")
        if au_count >= 3:
            score += 1
            details["tail_au"] = {"value": au_count, "passed": True}
        else:
            violations.append(f"Positions 15-19 have only {au_count} A/U (need >= 3)")
            details["tail_au"] = {"value": au_count, "passed": False}

        # 3. 无内部重复（Tm < 20°C）—— 简化为检查是否存在 >= 7bp 的反向重复
        # 实际 Tm 计算较复杂，这里用简化启发式
        has_repeat = self._check_internal_repeat(core)
        if not has_repeat:
            score += 1
            details["internal_repeat"] = {"passed": True}
        else:
            violations.append("Internal repeat detected (simplified check)")
            details["internal_repeat"] = {"passed": False}

        # 4. 19 位为 A
        if core[18] == "A":
            score += 1
            details["pos19_A"] = {"passed": True}
        else:
            violations.append("Position 19 is not A")
            details["pos19_A"] = {"passed": False}

        # 5. 3 位为 A
        if core[2] == "A":
            score += 1
            details["pos3_A"] = {"passed": True}
        else:
            violations.append("Position 3 is not A")
            details["pos3_A"] = {"passed": False}

        # 6. 10 位为 U
        if core[9] == "U":
            score += 1
            details["pos10_U"] = {"passed": True}
        else:
            violations.append("Position 10 is not U")
            details["pos10_U"] = {"passed": False}

        # 7. 19 位不为 G/C
        if core[18] not in "GC":
            score += 1
            details["pos19_not_GC"] = {"passed": True}
        else:
            # 注意：如果 19 位是 A，同时满足 #4 和 #7，但只加一次分
            # 这里逻辑需要和 #4 协调
            # 实际上 Reynolds 规则中 #4 和 #7 是独立的加分项
            # #4: 19 位是 A (+1), #7: 19 位不是 G/C (+1)
            # 如果 19 位是 A，同时满足两者，得 2 分
            # 如果 19 位是 U，满足 #7 不满足 #4，得 1 分
            # 如果 19 位是 G/C，两者都不满足，得 0 分
            pass  # 已在 #4 中处理，这里不做重复扣分

        # 修正 #7 逻辑
        if core[18] not in "GC":
            # 如果之前因为不是 A 没加到 #4 的分，但这里是 U，满足 #7
            if core[18] != "A":
                score += 1
                details["pos19_not_GC"] = {"passed": True}
        else:
            violations.append("Position 19 is G/C")
            details["pos19_not_GC"] = {"passed": False}

        # 8. 13 位不为 G
        if core[12] != "G":
            score += 1
            details["pos13_not_G"] = {"passed": True}
        else:
            violations.append("Position 13 is G")
            details["pos13_not_G"] = {"passed": False}

        passed = score >= self.CUTOFF

        return {
            "score": float(score),
            "passed": passed,
            "violations": violations,
            "details": details,
        }

    def _check_internal_repeat(self, sequence: str) -> bool:
        """
        简化检查：是否存在长度 >= 7 的反向互补子序列
        """
        complement = {"A": "U", "U": "A", "G": "C", "C": "G"}
        seq = sequence.upper()
        n = len(seq)

        for i in range(n):
            for length in range(7, min(12, n - i + 1)):
                sub = seq[i : i + length]
                rev_comp = "".join(complement.get(b, b) for b in reversed(sub))
                # 检查反向互补是否出现在序列中
                if rev_comp in seq and rev_comp != sub:
                    return True
        return False
