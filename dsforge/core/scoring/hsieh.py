"""
Hsieh 评分器 (2004)
- 22 个基因，138 条 siRNA
- 位置特征 + 热力学稳定性

简化实现：结合位置偏好和热力学特征
"""

from dsforge.core.scoring.base import ScoringRule, register_rule


@register_rule
class HsiehRule(ScoringRule):
    """
    Hsieh et al. (2004) — A new scheme for siRNA design
    """

    name = "hsieh"
    description = "Hsieh rules (2004): positional features + thermodynamic stability"

    # Hsieh 位置特征权重
    POSITION_WEIGHTS = {
        1:  {"A": -1, "U": -1, "G": +2, "C": +1},  # 5' 端 G 强偏好
        2:  {"A": -1, "U": -1, "G": +1, "C": +1},
        3:  {"A": +1, "U": -1, "G": +1, "C": -1},
        6:  {"A": +1, "U": +1, "G": -1, "C": -1},  # A/U 偏好
        10: {"A": +1, "U": +1, "G": 0,  "C": 0},
        13: {"A": +1, "U": +1, "G": -1, "C": -1},
        14: {"A": +1, "U": +1, "G": -1, "C": -1},
        18: {"A": +1, "U": +1, "G": -1, "C": -1},
        19: {"A": +2, "U": +2, "G": -2, "C": -2}, # 3' 端强 A/U 偏好
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

        # 位置特征评分
        pos_score = 0
        for pos, weights in self.POSITION_WEIGHTS.items():
            base = core[pos - 1]
            val = weights.get(base, 0)
            pos_score += val
            if val < 0:
                violations.append(f"Position {pos}: {base} has negative weight ({val})")

        details["positional_score"] = pos_score
        score += max(0, pos_score)

        # 热力学稳定性：5' 端 vs 3' 端不对称性
        # Hsieh 认为 antisense 5' 端更不稳定（更多 A/U）有利于 RISC 加载
        five_prime = core[:7]
        three_prime = core[12:19]
        five_au = five_prime.count("A") + five_prime.count("U")
        three_au = three_prime.count("A") + three_prime.count("U")
        asymmetry = five_au - three_au

        if asymmetry >= 2:
            score += 2
            details["asymmetry"] = {"value": asymmetry, "passed": True, "note": "5' more A/U than 3'"}
        else:
            violations.append(f"5'/3' asymmetry {asymmetry} (need >= 2)")
            details["asymmetry"] = {"value": asymmetry, "passed": False}

        # GC 含量 30-55%
        gc_pct = ((core.count("G") + core.count("C")) / 19) * 100
        if 30 <= gc_pct <= 55:
            score += 1
            details["gc_content"] = {"value": gc_pct, "passed": True}
        else:
            violations.append(f"GC content {gc_pct:.1f}% not in 30-55%")
            details["gc_content"] = {"value": gc_pct, "passed": False}

        # 归一化（理论最大约 15 分）
        normalized = min(100, (score / 10) * 100)
        passed = score >= 5

        return {
            "score": round(normalized, 2),
            "passed": passed,
            "violations": violations,
            "details": details,
        }
