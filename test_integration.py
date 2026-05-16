#!/usr/bin/env python3
"""
端到端集成测试
验证三种设计模式 + 导出 + 数据库完整流程
"""
import sys, os, tempfile, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dsforge.core.sequence import TranscriptomeIndex
from dsforge.controller.design_task import DesignTask, DesignConfig
from dsforge.controller.design_task_parallel import ParallelDesignTask
from dsforge.controller.exporter import ResultExporter
from dsforge.database.manager import DatabaseManager


def create_test_transcriptome():
    """创建测试转录组"""
    t = TranscriptomeIndex()
    t.sequences = {
        "gene_A": "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC",
        "gene_B": "GGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUU",
    }
    t._compute_stats()
    return t


def test_end_to_end_siRNA():
    print("=" * 60)
    print("E2E Test: siRNA Mode")
    print("=" * 60)

    t = create_test_transcriptome()
    db = DatabaseManager(":memory:")
    task = DesignTask(db_manager=db)

    config = DesignConfig(
        mode="siRNA",
        enabled_rules=["consensus", "reynolds", "ui_tei", "amarzguioui", "hsieh", "jagla"],
        gc_min=20, gc_max=80,
    )

    result = task.run(t, "gene_A", config)
    assert result["task_id"] == 1
    assert result["mode"] == "siRNA"
    assert len(result["results"]) > 0
    print(f"  ✓ {len(result['results'])} candidates, top score: {result['results'][0]['consensus_score']:.1f}")

    # 验证数据库
    task_record = db.get_task(result["task_id"])
    assert task_record["status"] == "completed"
    print(f"  ✓ Database record verified")

    # 导出测试
    exporter = ResultExporter()
    with tempfile.TemporaryDirectory() as tmpdir:
        exporter.export_csv(result["results"], f"{tmpdir}/results.csv")
        exporter.export_fasta(result["results"], f"{tmpdir}/results.fa")
        assert os.path.exists(f"{tmpdir}/results.csv")
        assert os.path.exists(f"{tmpdir}/results.fa")
        print(f"  ✓ Export verified")

    print("  ✓ siRNA E2E OK\n")


def test_end_to_end_long_dsRNA_parallel():
    print("=" * 60)
    print("E2E Test: Long dsRNA Mode (Parallel)")
    print("=" * 60)

    t = create_test_transcriptome()
    db = DatabaseManager(":memory:")
    task = ParallelDesignTask(db_manager=db, n_cores=2, batch_size=20)

    config = DesignConfig(
        mode="long_dsRNA",
        length_min=50, length_max=80,
        enabled_rules=["consensus", "reynolds", "ui_tei"],
        gc_min=20, gc_max=80,
    )

    result = task.run_parallel(t, "gene_A", config)
    assert result["mode"] == "long_dsRNA"
    assert len(result["results"]) > 0
    top = result["results"][0]
    print(f"  ✓ {len(result['results'])} regions")
    print(f"  ✓ Top region: {len(top['sequence'])} bp, score: {top['consensus_score']:.1f}")
    print(f"  ✓ Pool products: {top.get('pool', {}).get('num_products', 0)}")

    # 验证数据库持久化
    results = db.get_results(result["task_id"])
    assert len(results) > 0
    print(f"  ✓ Database: {len(results)} result rows")

    print("  ✓ Long dsRNA Parallel E2E OK\n")


def test_end_to_end_dsiRNA():
    print("=" * 60)
    print("E2E Test: DsiRNA Mode")
    print("=" * 60)

    t = create_test_transcriptome()
    db = DatabaseManager(":memory:")
    task = DesignTask(db_manager=db)

    config = DesignConfig(
        mode="DsiRNA",
        enabled_rules=["consensus", "reynolds"],
        gc_min=20, gc_max=80,
    )

    result = task.run(t, "gene_A", config)
    assert result["mode"] == "DsiRNA"
    assert len(result["results"]) > 0
    top = result["results"][0]
    print(f"  ✓ {len(result['results'])} candidates")
    print(f"  ✓ Top DsiRNA: {top['sequence']}, pool score: {top['consensus_score']:.1f}")
    print("  ✓ DsiRNA E2E OK\n")


def test_history_and_reload():
    print("=" * 60)
    print("E2E Test: History & Reload")
    print("=" * 60)

    t = create_test_transcriptome()
    db = DatabaseManager(":memory:")
    task = DesignTask(db_manager=db)

    # 运行两个任务
    task.run(t, "gene_A", DesignConfig(mode="siRNA", enabled_rules=["consensus"]))
    task.run(t, "gene_B", DesignConfig(mode="long_dsRNA", length_min=40, length_max=50, enabled_rules=["consensus"]))

    # 列出历史
    history = db.list_tasks()
    assert len(history) == 2
    print(f"  ✓ History: {len(history)} tasks")

    # 重载第一个任务的结果
    results = db.get_results(history[1]["id"])  # 最新的
    assert len(results) > 0
    print(f"  ✓ Reloaded {len(results)} results from task {history[1]['id']}")
    print("  ✓ History & Reload OK\n")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("dsRNA-Forge End-to-End Integration Test")
    print("=" * 60 + "\n")

    tests = [
        test_end_to_end_siRNA,
        test_end_to_end_long_dsRNA_parallel,
        test_end_to_end_dsiRNA,
        test_history_and_reload,
    ]

    passed = 0
    failed = 0
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            failed += 1
            print(f"❌ {test.__name__} FAILED: {e}")
            import traceback
            traceback.print_exc()

    print("=" * 60)
    if failed == 0:
        print(f"ALL {passed} E2E TESTS PASSED ✓")
    else:
        print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 60)

    sys.exit(failed)
