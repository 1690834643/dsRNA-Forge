"""
Jagla 评分器 (2005)
基于 601 条 siRNA 的 4 套决策树规则（Jagla1–Jagla4）
分别针对不同 GC 含量区间
"""

from dsforge.core.scoring.base import ScoringRule, register_rule


@register_rule
class JaglaRule(ScoringRule):
    """
    Jagla et al. (2005) — Sequence characteristics of functional siRNAs
    4 套决策树，按 GC 含量区间选择
    """

    name = "jagla"
    description = "Jagla rules (2005): 4 decision trees for different GC ranges"

    # Jagla 决策树规则（简化版，基于文献描述）
    # Jagla1: GC < 35%
    # Jagla2: 35% <= GC < 45%
    # Jagla3: 45% <= GC < 55%
    # Jagla4: GC >= 55%

    JAGLA_RULES = {
        1: {  # GC < 35%
            "gc_range": (0, 35),
            "pos1_pref": "GC",
            "pos19_pref": "AU",
            "min_tail_au": 4,  # 13-19 位 A/U 数
            "max_gc_run": 7,
        },
        2: {  # 35% <= GC < 45%
            "gc_range": (35, 45),
            "pos1_pref": "GC",
            "pos19_pref": "AU",
            "min_tail_au": 3,
            "max_gc_run": 8,
        },
        3: {  # 45% <= GC < 55%
            "gc_range": (45, 55),
            "pos1_pref": "GC",
            "pos19_pref": "AU",
            "min_tail_au": 3,
            "max_gc_run": 9,
        },
        4: {  # GC >= 55%
            "gc_range": (55, 100),
            "pos1_pref": "GC",
            "pos19_pref": "AU",
            "min_tail_au": 2,
            "max_gc_run": 10,
        },
    }

    def score(self, sequence: str) -> dict:
        if not self.validate_sequence(sequence):
            return {"score": 0.0, "passed": False, "violations": ["Invalid sequence"], "details": {}}

        seq = sequence.upper().replace("T", "U")
        if len(seq) < 19:
            return {"score": 0.0, "passed": False, "violations": ["Sequence too short (< 19nt)"], "details": {}}

        core = seq[:19]
        gc_pct = ((core.count("G") + core.count("C")) / 19) * 100

        # 选择适用的决策树
        jagla_id = self._select_jagla(gc_pct)
        rules = self.JAGLA_RULES[jagla_id]

        score = 0
        violations = []
        details = {
            "jagla_id": jagla_id,
            "gc_pct": gc_pct,
            "gc_range": rules["gc_range"],
        }

        # 1. 1 位偏好
        if core[0] in rules["pos1_pref"]:
            score += 1
            details["pos1"] = {"value": core[0], "passed": True}
        else:
            violations.append(f"Position 1 is {core[0]} (need {rules['pos1_pref']})")
            details["pos1"] = {"value": core[0], "passed": False}

        # 2. 19 位偏好
        if core[18] in rules["pos19_pref"]:
            score += 1
            details["pos19"] = {"value": core[18], "passed": True}
        else:
            violations.append(f"Position 19 is {core[18]} (need {rules['pos19_pref']})")
            details["pos19"] = {"value": core[18], "passed": False}

        # 3. 尾部 A/U 数
        tail = core[12:19]
        au_count = tail.count("A") + tail.count("U")
        if au_count >= rules["min_tail_au"]:
            score += 1
            details["tail_au"] = {"value": au_count, "passed": True, "min_required": rules["min_tail_au"]}
        else:
            violations.append(f"Positions 13-19 have {au_count} A/U (need >= {rules['min_tail_au']})")
            details["tail_au"] = {"value": au_count, "passed": False, "min_required": rules["min_tail_au"]}

        # 4. GC 连续长度
        max_gc_run = self._max_consecutive_gc(core)
        if max_gc_run <= rules["max_gc_run"]:
            score += 1
            details["max_gc_run"] = {"value": max_gc_run, "passed": True, "max_allowed": rules["max_gc_run"]}
        else:
            violations.append(f"GC run length {max_gc_run} (max allowed {rules['max_gc_run']})")
            details["max_gc_run"] = {"value": max_gc_run, "passed": False, "max_allowed": rules["max_gc_run"]}

        # 5. 整体 GC 范围检查
        if rules["gc_range"][0] <= gc_pct < rules["gc_range"][1]:
            score += 1
            details["gc_in_range"] = {"passed": True}
        else:
            # 这个条件通常已经由 jagla_id 选择保证
            pass

        # 归一化
        normalized = (score / 5) * 100
        passed = score >= 4

        return {
            "score": round(normalized, 2),
            "passed": passed,
            "violations": violations,
            "details": details,
        }

    def _select_jagla(self, gc_pct: float) -> int:
        """根据 GC 含量选择 Jagla 决策树"""
        for jagla_id, rules in self.JAGLA_RULES.items():
            low, high = rules["gc_range"]
            if low <= gc_pct < high:
                return jagla_id
        return 4  # 默认最高 GC 区间

    def _max_consecutive_gc(self, sequence: str) -> int:
        max_run = 0
        current = 0
        for base in sequence.upper():
            if base in "GC":
                current += 1
                max_run = max(max_run, current)
            else:
                current = 0
        return max_run
