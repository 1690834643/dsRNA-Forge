"""
共识评分器（Consensus Scorer）
基于 PMC5357899 (2017) 和常用 siRNA 规则重构的启发式位置编码表

对每个位置（1-19）的每个碱基（A/C/G/U），根据多数规则的共识赋予编码：
  +1 : 多数规则认为该碱基在此位置偏好高效 siRNA
  -1 : 多数规则认为该碱基在此位置偏好低效 siRNA
   0 : 无共识
"""

from dsforge.core.scoring.base import ScoringRule, register_rule


# 共识编码表（位置 1-19，碱基 A/C/G/U）
# 来源：PMC5357899 (2017) + 调研报告中的规则汇总。
# 注：这是离线筛选用的重构启发式，不应宣传为原文完整模型。
CONSENSUS_TABLE = {
    # 位置: {A: 编码, C: 编码, G: 编码, U: 编码}
    1:  {"A": -1, "C": +1, "G": +1, "U": -1},   # Ui-Tei, Amarzguioui, Jagla1-3 偏好 C/G
    2:  {"A": -1, "C":  0, "G": +1, "U": -1},   # Amarzguioui, Svetlana 偏好 G
    3:  {"A": +1, "C": -1, "G": +1, "U": -1},   # Reynolds, Svetlana, Jiang 偏好 A/G
    4:  {"A":  0, "C":  0, "G":  0, "U":  0},   # 无明确共识
    5:  {"A":  0, "C":  0, "G":  0, "U":  0},   # 无明确共识
    6:  {"A": +1, "C": -1, "G": -1, "U": +1},   # Amarzguioui, Hsieh, Takasaki 偏好 A/U
    7:  {"A":  0, "C":  0, "G":  0, "U":  0},   # 无明确共识
    8:  {"A":  0, "C":  0, "G":  0, "U":  0},   # 无明确共识
    9:  {"A":  0, "C":  0, "G":  0, "U":  0},   # 无明确共识
    10: {"A": +1, "C": +1, "G": +1, "U": +1},   # Reynolds, Svetlana 等大多数位置偏好
    11: {"A":  0, "C":  0, "G":  0, "U":  0},   # 无明确共识
    12: {"A":  0, "C":  0, "G":  0, "U":  0},   # 无明确共识
    13: {"A": +1, "C": -1, "G": -1, "U": +1},   # Reynolds: 13 位不为 G
    14: {"A":  0, "C":  0, "G":  0, "U":  0},   # 无明确共识
    15: {"A": +1, "C": -1, "G": -1, "U": +1},   # 尾部偏好 A/U
    16: {"A": +1, "C": -1, "G": -1, "U": +1},   # 尾部偏好 A/U
    17: {"A": +1, "C": -1, "G": -1, "U": +1},   # 尾部偏好 A/U
    18: {"A": +1, "C": -1, "G": -1, "U": +1},   # 尾部偏好 A/U
    19: {"A": +1, "C": -1, "G": -1, "U": +1},   # Ui-Tei, Amarzguioui, Hsieh 等一致偏好 A/U
}


@register_rule
class ConsensusRule(ScoringRule):
    """
    共识评分器 — 默认引擎
    基于 PMC5357899 (2017) 和常用规则重构的启发式位置编码表
    """

    name = "consensus"
    description = "Heuristic consensus scorer inspired by PMC5357899 (2017)"

    def score(self, sequence: str) -> dict:
        """
        对 21nt siRNA antisense 序列进行共识评分
        """
        if not self.validate_sequence(sequence):
            return {"score": 0.0, "passed": False, "violations": ["Invalid sequence"], "details": {}}

        seq = sequence.upper().replace("T", "U")
        if len(seq) < 19:
            return {"score": 0.0, "passed": False, "violations": ["Sequence too short (< 19nt)"], "details": {}}

        core = seq[:19]
        total_score = 0
        position_scores = []
        violations = []

        for i, base in enumerate(core, start=1):
            code = CONSENSUS_TABLE.get(i, {}).get(base, 0)
            total_score += code
            position_scores.append({
                "position": i,
                "base": base,
                "code": code,
            })
            if code == -1:
                violations.append(f"Position {i}: {base} is consensus-unfavorable")

        # 归一化到 0-100 范围便于比较
        # 理论最大 +19，最小 -19
        normalized_score = ((total_score + 19) / 38) * 100

        # 共识评分器的 pass 逻辑：正分即通过
        passed = total_score > 0

        return {
            "score": round(normalized_score, 2),
            "passed": passed,
            "violations": violations,
            "details": {
                "raw_score": total_score,
                "max_possible": 19,
                "min_possible": -19,
                "position_scores": position_scores,
            },
        }
