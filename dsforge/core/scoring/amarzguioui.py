"""
Amarzguioui 评分器 (2004)
- 位置特异性碱基偏好
- 无内部重复
- 6 位为 A 偏好高效
"""

from dsforge.core.scoring.base import ScoringRule, register_rule


@register_rule
class AmarzguiouiRule(ScoringRule):
    """
    Amarzguioui et al. (2004) — Tolerance for mutations and chemical modifications in a siRNA
    """

    name = "amarzguioui"
    description = "Amarzguioui rules (2004): positional preferences + internal repeat"

    # 位置偏好编码（基于文献中的位置特异性偏好）
    POSITION_PREFERENCES = {
        1:  {"A": 0, "C": +1, "G": +1, "U": 0},   # 5' 端偏好 G/C
        6:  {"A": +1, "C": 0, "G": 0, "U": 0},    # 6 位 A 偏好
        10: {"A": 0, "C": 0, "G": 0, "U": +1},    # 10 位 U 偏好
        13: {"A": +1, "C": -1, "G": -1, "U": +1}, # 13 位 A/U 偏好
        19: {"A": +1, "C": -1, "G": -1, "U": +1}, # 3' 端 A/U 偏好
    }

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

        # 位置特异性评分
        pos_score = 0
        for pos, prefs in self.POSITION_PREFERENCES.items():
            base = core[pos - 1]
            val = prefs.get(base, 0)
            pos_score += val
            if val < 0:
                violations.append(f"Position {pos}: {base} is unfavorable")

        details["positional_score"] = pos_score
        score += max(0, pos_score)  # 只加正分

        # 无内部重复（简化：检查 >= 7bp 反向互补）
        has_repeat = self._check_internal_repeat(core)
        if not has_repeat:
            score += 2
            details["internal_repeat"] = {"passed": True}
        else:
            violations.append("Internal repeat detected")
            details["internal_repeat"] = {"passed": False}

        # 6 位为 A
        if len(core) >= 6 and core[5] == "A":
            score += 1
            details["pos6_A"] = {"passed": True}
        else:
            violations.append("Position 6 is not A")
            details["pos6_A"] = {"passed": False}

        # 归一化（理论最大约 8 分）
        normalized = min(100, (score / 6) * 100)
        passed = score >= 3

        return {
            "score": round(normalized, 2),
            "passed": passed,
            "violations": violations,
            "details": details,
        }

    def _check_internal_repeat(self, sequence: str) -> bool:
        """检查是否存在长度 >= 7 的反向互补子序列"""
        complement = {"A": "U", "U": "A", "G": "C", "C": "G"}
        seq = sequence.upper()
        n = len(seq)
        for i in range(n):
            for length in range(7, min(12, n - i + 1)):
                sub = seq[i : i + length]
                rev_comp = "".join(complement.get(b, b) for b in reversed(sub))
                if rev_comp in seq and rev_comp != sub:
                    return True
        return False
