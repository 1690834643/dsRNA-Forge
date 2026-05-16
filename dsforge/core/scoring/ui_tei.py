"""
Ui-Tei 评分器 (2004)
基于 72 条 siRNA（靶向 6 个基因）的 4 条标准
"""

from dsforge.core.scoring.base import ScoringRule, register_rule


@register_rule
class UiTeiRule(ScoringRule):
    """
    Ui-Tei et al. (2004) — Guidelines for the selection of highly effective siRNA sequences
    """

    name = "ui_tei"
    description = "Ui-Tei rules (2004): 4 criteria for antisense strand"

    def score(self, sequence: str) -> dict:
        if not self.validate_sequence(sequence):
            return {"score": 0.0, "passed": False, "violations": ["Invalid sequence"], "details": {}}

        seq = sequence.upper().replace("T", "U")
        if len(seq) < 19:
            return {"score": 0.0, "passed": False, "violations": ["Sequence too short (< 19nt)"], "details": {}}

        core = seq[:19]
        score = 0
        violations = []
        details = {}

        # 1. 19 位（antisense 3' 端）为 A 或 U
        if core[18] in "AU":
            score += 1
            details["pos19_AU"] = {"value": core[18], "passed": True}
        else:
            violations.append(f"Position 19 is {core[18]} (need A/U)")
            details["pos19_AU"] = {"value": core[18], "passed": False}

        # 2. 1 位（antisense 5' 端）为 G 或 C
        if core[0] in "GC":
            score += 1
            details["pos1_GC"] = {"value": core[0], "passed": True}
        else:
            violations.append(f"Position 1 is {core[0]} (need G/C)")
            details["pos1_GC"] = {"value": core[0], "passed": False}

        # 3. 13-19 位至少 5 个 A/U
        tail = core[12:19]  # 位置 13-19
        au_count = tail.count("A") + tail.count("U")
        if au_count >= 5:
            score += 1
            details["tail13_19_au"] = {"value": au_count, "passed": True}
        else:
            violations.append(f"Positions 13-19 have {au_count} A/U (need >= 5)")
            details["tail13_19_au"] = {"value": au_count, "passed": False}

        # 4. GC 连续不超过 9 nt
        max_gc_run = self._max_consecutive_gc(core)
        if max_gc_run <= 9:
            score += 1
            details["max_gc_run"] = {"value": max_gc_run, "passed": True}
        else:
            violations.append(f"GC run length {max_gc_run} (max allowed 9)")
            details["max_gc_run"] = {"value": max_gc_run, "passed": False}

        # 归一化：4 条标准，每条 25 分
        normalized = (score / 4) * 100
        passed = score >= 3  # 至少满足 3/4

        return {
            "score": round(normalized, 2),
            "passed": passed,
            "violations": violations,
            "details": details,
        }

    def _max_consecutive_gc(self, sequence: str) -> int:
        """计算最长连续 GC 长度"""
        max_run = 0
        current = 0
        for base in sequence.upper():
            if base in "GC":
                current += 1
                max_run = max(max_run, current)
            else:
                current = 0
        return max_run
