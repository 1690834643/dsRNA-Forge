#!/usr/bin/env python3
"""
CLI 测试脚本 — 验证核心引擎（含 ViennaRNA）
"""

import sys
import os

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, PROJECT_ROOT)

from dsforge.core.sequence import TranscriptomeIndex, generate_candidates, gc_content
from dsforge.core.scoring.base import evaluate_all_rules, list_rules
from dsforge.core.scoring.reynolds import ReynoldsRule
from dsforge.core.scoring.consensus import ConsensusRule
from dsforge.core.scoring.ui_tei import UiTeiRule
from dsforge.core.scoring.amarzguioui import AmarzguiouiRule
from dsforge.core.scoring.hsieh import HsiehRule
from dsforge.core.scoring.jagla import JaglaRule
from dsforge.core.dicer import predict_dicer_products, calculate_pool_score
from dsforge.core.offtarget import OffTargetScreener
from dsforge.core.thermodynamics import ThermodynamicsCalculator
from dsforge.database.manager import DatabaseManager
from dsforge.controller.design_task import DesignTask, DesignConfig


def test_sequence_utils():
    print("=" * 60)
    print("Test 1: Sequence Utilities")
    print("=" * 60)

    seq = "AUGCAUGCAUGCAUGCAUGC"
    print(f"Sequence: {seq}")
    print(f"GC content: {gc_content(seq):.1f}%")

    candidates = list(generate_candidates(seq, "siRNA", 5, 7, gc_min=30, gc_max=70))
    print(f"Generated {len(candidates)} candidates (len 5-7)")
    for c in candidates[:3]:
        print(f"  {c['sequence']} @ {c['start']}-{c['end']} GC={c['gc']:.1f}%")

    print("✓ Sequence utilities OK\n")


def test_scoring():
    print("=" * 60)
    print("Test 2: Scoring Engine (All 6 Rules)")
    print("=" * 60)

    print(f"Registered rules: {list_rules()}")
    assert len(list_rules()) == 6, f"Expected 6 rules, got {len(list_rules())}"

    test_seq = "AAGGCUAUGUAGAUUUAUGCC"

    # 单独测试每个规则
    rules = {
        "reynolds": ReynoldsRule(),
        "consensus": ConsensusRule(),
        "ui_tei": UiTeiRule(),
        "amarzguioui": AmarzguiouiRule(),
        "hsieh": HsiehRule(),
        "jagla": JaglaRule(),
    }

    print(f"\nTest sequence: {test_seq}")
    for name, rule in rules.items():
        result = rule.score(test_seq)
        status = "PASS" if result['passed'] else "FAIL"
        print(f"  {name:15s}: {result['score']:6.1f} [{status}] violations={len(result['violations'])}")

    # 多规则评估
    results = evaluate_all_rules(test_seq, list(rules.keys()))
    print(f"\nMulti-rule evaluation:")
    for rule, res in results.items():
        print(f"  {rule:15s}: {res['score']:6.1f} (passed={res['passed']})")

    print("✓ All 6 scoring rules OK\n")


def test_dicer():
    print("=" * 60)
    print("Test 3: Dicer Cleavage Prediction")
    print("=" * 60)

    dsRNA = "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC"
    products = predict_dicer_products(dsRNA, cut_length=21)
    print(f"dsRNA length: {len(dsRNA)}")
    print(f"Predicted Dicer products: {len(products)}")

    for p in products[:3]:
        print(f"  {p['sequence']} @ pos {p['position']} len={p['length']} +3G={p['plus3g_bonus']}")

    pool = calculate_pool_score(products, evaluate_all_rules, ["consensus"])
    print(f"\nPool score: {pool['pool_score']:.2f}")
    print(f"  Products: {pool['num_products']}, High quality: {pool['high_quality_count']}")

    print("✓ Dicer prediction OK\n")


def test_database():
    print("=" * 60)
    print("Test 4: Database")
    print("=" * 60)

    db = DatabaseManager(":memory:")

    task_id = db.create_task(
        mode="siRNA",
        target_seq_id="test_gene",
        target_seq="AUGCAUGCAUGC",
        params={"test": True},
    )
    print(f"Created task: {task_id}")

    db.update_task_status(task_id, "running")
    db.update_task_status(task_id, "completed")

    task = db.get_task(task_id)
    print(f"Task status: {task['status']}")

    result_id = db.add_result(task_id, 1, "AUGCAUGCAUGC", 0, 12, 85.5, 1)
    print(f"Added result: {result_id}")

    db.add_rule_score(result_id, "consensus", 85.5, True, [])
    db.add_thermodynamics(result_id, -15.2, 0, 0, -18.5)

    results = db.get_results(task_id)
    print(f"Retrieved {len(results)} result(s)")

    print("✓ Database OK\n")


def test_design_task_siRNA():
    print("=" * 60)
    print("Test 5: Design Task (siRNA mode)")
    print("=" * 60)

    transcriptome = TranscriptomeIndex()
    transcriptome.sequences = {
        "test_gene_1": "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC",
        "test_gene_2": "GGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUU",
    }
    transcriptome._compute_stats()
    print(f"Transcriptome: {transcriptome.get_stats()}")

    config = DesignConfig(
        mode="siRNA",
        enabled_rules=["consensus", "reynolds", "ui_tei"],
        gc_min=20,
        gc_max=80,
    )

    task = DesignTask(db_manager=DatabaseManager(":memory:"))

    def progress_cb(step, pct):
        if int(pct) % 20 == 0:
            print(f"  [{pct:5.1f}%] {step}")

    result = task.run(
        transcriptome=transcriptome,
        target_seq_id="test_gene_1",
        config=config,
        progress_callback=progress_cb,
    )

    print(f"\nTask completed: {result['task_id']}")
    print(f"Total candidates: {result['summary']['total_candidates']}")
    print(f"Passed: {result['summary']['passed_candidates']}")

    if result['results']:
        top = result['results'][0]
        print(f"\nTop candidate:")
        print(f"  Rank: {top['rank']}")
        print(f"  Sequence: {top['sequence']}")
        print(f"  Score: {top['consensus_score']:.2f}")
        print(f"  Passed: {top['passed']}")

    print("✓ siRNA mode OK\n")


def test_design_task_long_dsRNA():
    print("=" * 60)
    print("Test 6: Design Task (Long dsRNA mode)")
    print("=" * 60)

    # 使用更长的测试序列
    long_seq = "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC"

    transcriptome = TranscriptomeIndex()
    transcriptome.sequences = {
        "test_gene_long": long_seq,
    }
    transcriptome._compute_stats()
    print(f"Target length: {len(long_seq)} nt")

    config = DesignConfig(
        mode="long_dsRNA",
        length_min=50,
        length_max=80,
        enabled_rules=["consensus", "reynolds", "ui_tei", "jagla"],
        gc_min=20,
        gc_max=80,
    )

    task = DesignTask(db_manager=DatabaseManager(":memory:"))

    def progress_cb(step, pct):
        if int(pct) % 20 == 0 or pct >= 99:
            print(f"  [{pct:5.1f}%] {step}")

    result = task.run(
        transcriptome=transcriptome,
        target_seq_id="test_gene_long",
        config=config,
        progress_callback=progress_cb,
    )

    print(f"\nTask completed: {result['task_id']}")
    print(f"Total regions: {result['summary']['total_candidates']}")
    print(f"Passed: {result['summary']['passed_candidates']}")

    if result['results']:
        top = result['results'][0]
        print(f"\nTop region:")
        print(f"  Rank: {top['rank']}")
        print(f"  Length: {len(top['sequence'])} bp")
        print(f"  Position: {top['position']}")
        print(f"  Pool Score: {top['consensus_score']:.2f}")
        print(f"  Passed: {top['passed']}")
        pool = top.get('pool', {})
        print(f"  Dicer products: {pool.get('num_products', 0)}")
        print(f"  High-quality products: {pool.get('high_quality_count', 0)}")

    print("✓ Long dsRNA mode OK\n")


def test_design_task_dsiRNA():
    print("=" * 60)
    print("Test 7: Design Task (DsiRNA mode)")
    print("=" * 60)

    transcriptome = TranscriptomeIndex()
    transcriptome.sequences = {
        "test_gene_1": "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC",
    }
    transcriptome._compute_stats()

    config = DesignConfig(
        mode="DsiRNA",
        enabled_rules=["consensus", "reynolds"],
        gc_min=20,
        gc_max=80,
    )

    task = DesignTask(db_manager=DatabaseManager(":memory:"))

    def progress_cb(step, pct):
        if int(pct) % 20 == 0 or pct >= 99:
            print(f"  [{pct:5.1f}%] {step}")

    result = task.run(
        transcriptome=transcriptome,
        target_seq_id="test_gene_1",
        config=config,
        progress_callback=progress_cb,
    )

    print(f"\nTask completed: {result['task_id']}")
    print(f"Total candidates: {result['summary']['total_candidates']}")
    print(f"Passed: {result['summary']['passed_candidates']}")

    if result['results']:
        top = result['results'][0]
        print(f"\nTop DsiRNA:")
        print(f"  Rank: {top['rank']}")
        print(f"  Sequence: {top['sequence']}")
        print(f"  Pool Score: {top['consensus_score']:.2f}")
        print(f"  Passed: {top['passed']}")

    print("✓ DsiRNA mode OK\n")


def test_thermodynamics():
    print("=" * 60)
    print("Test 8: ViennaRNA Thermodynamics")
    print("=" * 60)

    calc = ThermodynamicsCalculator()
    print(f"ViennaRNA available: {calc.available}")

    if calc.available:
        seq1, seq2 = 'GGGAAACCC', 'GGGUUUCCC'

        d = calc.rnaduplex(seq1, seq2)
        print(f"RNAduplex: {d}")

        c = calc.rnacofold(seq1, seq2)
        print(f"RNAcofold: {c}")

        u = calc.rnaup(seq1, seq2)
        print(f"RNAup (fallback): {u}")

        sd = calc.calculate_seed_dg('GAAACCC', 'AAAGGGUUUCCC')
        print(f"Seed DG: {sd}")
    else:
        print("  (skipped — ViennaRNA not installed)")

    print("✓ Thermodynamics OK\n")


def test_exporter():
    print("=" * 60)
    print("Test 9: Result Exporter")
    print("=" * 60)

    from dsforge.controller.exporter import ResultExporter

    exporter = ResultExporter()
    results = [
        {"rank": 1, "sequence": "AUGCAUGCAUGC", "position": "0-12", "consensus_score": 85.5, "passed": True},
        {"rank": 2, "sequence": "GCAUGCAUGCAU", "position": "1-13", "consensus_score": 72.3, "passed": True},
    ]

    exporter.export_csv(results, "/tmp/test_export.csv")
    print("CSV export: /tmp/test_export.csv")

    exporter.export_fasta(results, "/tmp/test_export.fa")
    print("FASTA export: /tmp/test_export.fa")

    print("✓ Exporter OK\n")


def main():
    print("\n" + "=" * 60)
    print("dsRNA-Forge CLI Test Suite — Phase 3 Verification")
    print("=" * 60 + "\n")

    tests = [
        test_sequence_utils,
        test_scoring,
        test_dicer,
        test_database,
        test_design_task_siRNA,
        test_design_task_long_dsRNA,
        test_design_task_dsiRNA,
        test_thermodynamics,
        test_exporter,
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            failed += 1
            print(f"\n❌ {test.__name__} FAILED: {e}")
            import traceback
            traceback.print_exc()
            print()

    print("=" * 60)
    if failed == 0:
        print(f"ALL {passed} TESTS PASSED ✓")
    else:
        print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 60)

    return failed


if __name__ == "__main__":
    sys.exit(main())
