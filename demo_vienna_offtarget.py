#!/usr/bin/env python3
"""
ViennaRNA 脱靶预测 Demo
展示热力学计算在实际转录组上的效果
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dsforge.core.sequence import TranscriptomeIndex
from dsforge.core.offtarget import OffTargetScreener
from dsforge.core.scoring.base import evaluate_all_rules

print("=" * 70)
print("ViennaRNA 脱靶预测 Demo")
print("=" * 70)

# 创建有脱靶关系的测试转录组
transcriptome = TranscriptomeIndex()
transcriptome.sequences = {
    # 目标基因
    "GeneA_target": "AUGCGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUU",
    # 脱靶基因：与 GeneA 有 16bp+ 连续匹配
    "GeneB_offtarget_16bp": "AAAAAAAUGCGGAAACCCUUUGGGAAACCCUUUCCCGGGUUUAAACCCGGGUUUAAACCCGGGUUU",
    # 脱靶基因：与 GeneA 种子区完美互补（用于 ViennaRNA 测试）
    "GeneC_seed_complement": "UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU",
}
transcriptome._compute_stats()

# 设计一条 siRNA，其种子区为 "GGAAACC"
test_sirna = "AUGCGGAAACCCUUUGGGAAA"

print(f"\n测试 siRNA: {test_sirna}")
print(f"种子区 (nt 2-8): {test_sirna[1:8]}")
print(f"转录组: {len(transcriptome.sequences)} 个基因")

# === 评分 ===
print("\n" + "-" * 70)
print("多规则评分:")
print("-" * 70)
results = evaluate_all_rules(test_sirna, ["consensus", "reynolds", "ui_tei"])
for rule_name, result in results.items():
    print(f"  {rule_name:12s}: score={result['score']:6.1f} passed={result['passed']}")

# === 脱靶筛查（带 ViennaRNA 热力学） ===
print("\n" + "-" * 70)
print("脱靶筛查（含 ViennaRNA RNAduplex 热力学）:")
print("-" * 70)

screener = OffTargetScreener(transcriptome)

start = time.time()
off_target = screener.screen_sequence(
    test_sirna,
    level_1_16bp=True,
    level_2_20bp=True,
    seed_7nt=True,
    use_vienna=True,
    seed_dg_threshold=-7.0,
)
elapsed = time.time() - start

print(f"\n  耗时: {elapsed*1000:.1f} ms")
print(f"  ViennaRNA 可用: {off_target['thermodynamics']['vienna_available']}")
print(f"  风险等级: {off_target['risk_level']}")
print(f"  是否通过: {off_target['passed']}")

# 序列比对脱靶
print(f"\n  序列比对脱靶:")
print(f"    16bp 连续匹配命中: {off_target['summary']['level_1_16bp_hits']}")
print(f"    20bp 连续匹配命中: {off_target['summary']['level_2_20bp_hits']}")
print(f"    7nt 种子区命中: {off_target['summary']['seed_7nt_hits']}")

# 热力学脱靶
thermo = off_target['thermodynamics']
print(f"\n  ViennaRNA 热力学脱靶:")
print(f"    最小 ΔG: {thermo['min_dg']} kcal/mol")
print(f"    热力学脱靶命中数: {len(thermo['seed_hits'])}")

if thermo['seed_hits']:
    print(f"\n    热力学脱靶详情:")
    for hit in thermo['seed_hits']:
        print(f"      → {hit['target_id']}: ΔG={hit['dg']:.2f} @ pos {hit['position']}")
        print(f"        seed: {hit['seed']} | target: {hit['target_seq']}")
        print(f"        structure: {hit.get('structure', 'N/A')}")
else:
    print(f"\n    无热力学脱靶（种子区结合能均高于阈值 {off_target.get('seed_dg_threshold', -7.0)} kcal/mol）")

# === Pool 脱靶筛查 ===
print("\n" + "-" * 70)
print("长 dsRNA Pool 脱靶筛查:")
print("-" * 70)

from dsforge.core.dicer import predict_dicer_products
products = predict_dicer_products(transcriptome.sequences["GeneA_target"][:60])
print(f"  Dicer 产物数: {len(products)}")

pool_offtarget = screener.screen_pool(products)
print(f"  Pool 风险: {pool_offtarget['pool_risk']}")
print(f"  高风险产物: {pool_offtarget['high_risk_products']}/{len(products)}")

# === 索引缓存测试 ===
print("\n" + "-" * 70)
print("转录组索引缓存测试:")
print("-" * 70)

# 先清除缓存
TranscriptomeIndex.clear_cache()

# 第一次加载（无缓存）
t1 = TranscriptomeIndex()
start = time.time()
t1.load_fasta("demo_data/test_transcriptome.fa", use_cache=True)
t1_elapsed = time.time() - start
print(f"  首次加载: {t1_elapsed*1000:.1f} ms (从 FASTA 解析)")

# 第二次加载（有缓存）
t2 = TranscriptomeIndex()
start = time.time()
t2.load_fasta("demo_data/test_transcriptome.fa", use_cache=True)
t2_elapsed = time.time() - start
print(f"  二次加载: {t2_elapsed*1000:.1f} ms (从缓存读取)")

if t1_elapsed > 0:
    speedup = t1_elapsed / t2_elapsed
    print(f"  加速比: {speedup:.1f}x")

cache_info = TranscriptomeIndex.get_cache_info()
print(f"  缓存文件: {len(cache_info)} 个")
for ci in cache_info:
    print(f"    {ci['file']}: {ci['size_mb']} MB")

print("\n" + "=" * 70)
print("Demo 完成 ✓")
print("=" * 70)
