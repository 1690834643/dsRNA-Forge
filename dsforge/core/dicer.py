"""
Dicer 切割位点预测器

Dicer 切割偏好：
- PAZ 结构域识别 dsRNA 3' 末端 2nt overhang
- 切割位点距 3' 末端约 21-22 nt
- +3G 偏好：第 3 位为 G 时 Dicer 结合更稳定
- 产物长度分布：主峰 21-22nt，拖尾 19-23nt

当前实现是用于候选相对排序的简化模拟，不是定量 Dicer 生化模型。
"""

from typing import List, Dict

# 默认切割参数
DEFAULT_CUT_LENGTH = 21
DEFAULT_CUT_VARIATION = 1  # ±1 nt (19-23nt)


def predict_dicer_products(
    dsRNA_seq: str,
    cut_length: int = DEFAULT_CUT_LENGTH,
    variation: int = DEFAULT_CUT_VARIATION,
    strand: str = "antisense",
) -> List[Dict]:
    """
    预测长 dsRNA 的 Dicer 切割产物

    Args:
        dsRNA_seq: dsRNA 的 antisense 链序列（5'->3'）
        cut_length: 默认切割长度（21nt）
        variation: 长度变异范围（±1nt）
        strand: "antisense" 或 "sense"

    Returns:
        产物列表，每项包含：
        {
            'sequence': str,        # 产物序列（antisense）
            'position': int,        # 在 dsRNA 中的起始位置
            'length': int,          # 产物长度
            'cut_site': int,        # Dicer 切割位点（距 5' 端）
            'plus3g_bonus': float,  # +3G 权重加分
        }
    """
    seq = dsRNA_seq.upper().replace("T", "U")
    products = []

    # 从 5' 端开始，以 cut_length 为步长扫描
    # 注意：Dicer 从 dsRNA 末端开始切割，产生 21-22nt 产物
    # 简化模型：固定步长 sliding window
    pos = 0
    while pos + 19 <= len(seq):
        # 产物长度有变异：19-23nt
        for length in range(cut_length - variation, cut_length + variation + 1):
            if pos + length > len(seq):
                continue

            product_seq = seq[pos : pos + length]

            # +3G 偏好权重
            plus3g_bonus = 0.0
            if len(product_seq) >= 3 and product_seq[2] == "G":
                plus3g_bonus = 0.5  # 权重加分

            products.append({
                "sequence": product_seq,
                "position": pos,
                "length": length,
                "cut_site": pos + length,
                "plus3g_bonus": plus3g_bonus,
            })

        pos += cut_length

    return products


def calculate_pool_score(
    products: List[Dict],
    scorer,
    enabled_rules: List[str],
    weight_by_plus3g: bool = True,
) -> Dict:
    """
    计算 siRNA Pool 的整体评分

    Args:
        products: Dicer 产物列表
        scorer: 评分函数（如 evaluate_all_rules）
        enabled_rules: 启用的规则列表
        weight_by_plus3g: 是否按 +3G 偏好加权

    Returns:
        {
            'pool_score': float,      # Pool 整体得分 (0-100)
            'num_products': int,      # 产物数量
            'avg_score': float,       # 平均得分
            'weighted_avg': float,    # 加权平均得分
            'high_quality_count': int,# 高质量产物数（passed=True）
            'low_quality_count': int, # 低质量产物数
            'product_details': List,  # 每个产物的详细评分
        }
    """
    if not products:
        return {
            "pool_score": 0.0,
            "num_products": 0,
            "avg_score": 0.0,
            "weighted_avg": 0.0,
            "high_quality_count": 0,
            "low_quality_count": 0,
            "product_details": [],
        }

    product_scores = []
    total_score = 0.0
    total_weight = 0.0
    high_quality = 0
    low_quality = 0

    for product in products:
        seq = product["sequence"]
        rules_result = scorer(seq, enabled_rules)

        # 取共识评分器的结果作为该产物的代表分
        consensus = rules_result.get("consensus", {})
        score = consensus.get("score", 0.0)
        passed = consensus.get("passed", False)

        # 权重
        weight = 1.0
        if weight_by_plus3g:
            weight += product.get("plus3g_bonus", 0.0)

        total_score += score * weight
        total_weight += weight

        if passed:
            high_quality += 1
        else:
            low_quality += 1

        product_scores.append({
            **product,
            "consensus_score": score,
            "passed": passed,
            "rule_details": rules_result,
        })

    avg_score = total_score / total_weight if total_weight > 0 else 0.0

    # Pool 整体得分 = 加权平均分 × (高质量比例)
    quality_ratio = high_quality / len(products) if products else 0
    pool_score = avg_score * quality_ratio

    return {
        "pool_score": round(pool_score, 2),
        "num_products": len(products),
        "avg_score": round(avg_score, 2),
        "weighted_avg": round(avg_score, 2),
        "high_quality_count": high_quality,
        "low_quality_count": low_quality,
        "product_details": product_scores,
    }
