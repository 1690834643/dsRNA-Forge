"""
多规则评分引擎基类与注册表
策略模式实现
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Any
import importlib


class ScoringRule(ABC):
    """
    评分规则基类
    所有具体规则（Reynolds, Ui-Tei 等）必须继承此类
    """

    name: str = ""
    description: str = ""

    @abstractmethod
    def score(self, sequence: str) -> Dict[str, Any]:
        """
        对一条 siRNA 序列（通常为 21nt antisense）进行评分

        Args:
            sequence: RNA 序列（5'->3'，antisense 链）

        Returns:
            {
                'score': float,          # 总分或归一化得分
                'passed': bool,          # 是否通过该规则的 cutoff
                'violations': List[str], # 未满足的具体条件（用于 GUI 显示）
                'details': Dict,         # 规则特有的详细指标
            }
        """
        pass

    def validate_sequence(self, sequence: str) -> bool:
        """验证序列是否合法（允许标准碱基 AUGC 及 IUPAC 歧义碱基）"""
        valid_bases = set("AUGCNRYMKSWHBVD")
        return set(sequence.upper()).issubset(valid_bases)


# 全局规则注册表
RULE_REGISTRY: Dict[str, type] = {}
_BUILTIN_RULES_LOADED = False
_BUILTIN_RULE_MODULES = [
    "dsforge.core.scoring.consensus",
    "dsforge.core.scoring.reynolds",
    "dsforge.core.scoring.ui_tei",
    "dsforge.core.scoring.amarzguioui",
    "dsforge.core.scoring.hsieh",
    "dsforge.core.scoring.jagla",
]


def register_rule(rule_class: type):
    """装饰器：将规则类注册到全局注册表"""
    if not issubclass(rule_class, ScoringRule):
        raise TypeError(f"{rule_class.__name__} must inherit from ScoringRule")
    RULE_REGISTRY[rule_class.name] = rule_class
    return rule_class


def ensure_builtin_rules_loaded():
    """Import bundled scoring rules so module-level decorators populate registry."""
    global _BUILTIN_RULES_LOADED
    if _BUILTIN_RULES_LOADED:
        return
    for module_name in _BUILTIN_RULE_MODULES:
        importlib.import_module(module_name)
    _BUILTIN_RULES_LOADED = True


def get_rule(name: str) -> ScoringRule:
    """通过名称获取规则实例"""
    ensure_builtin_rules_loaded()
    if name not in RULE_REGISTRY:
        raise KeyError(f"Rule '{name}' not found in registry. Available: {list(RULE_REGISTRY.keys())}")
    return RULE_REGISTRY[name]()


def list_rules() -> List[str]:
    """列出所有已注册的规则名称"""
    ensure_builtin_rules_loaded()
    return list(RULE_REGISTRY.keys())


def evaluate_all_rules(sequence: str, enabled_rules: List[str]) -> Dict[str, Dict[str, Any]]:
    """
    对一条序列运行所有启用的规则

    Args:
        sequence: RNA 序列
        enabled_rules: 启用的规则名称列表

    Returns:
        {rule_name: score_result, ...}
    """
    results = {}
    for rule_name in enabled_rules:
        try:
            rule = get_rule(rule_name)
            results[rule_name] = rule.score(sequence)
        except Exception as e:
            results[rule_name] = {
                "score": 0.0,
                "passed": False,
                "violations": [f"Evaluation error: {e}"],
                "details": {},
            }
    return results
