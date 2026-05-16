#!/usr/bin/env python3
"""多进程集成测试"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dsforge.core.sequence import TranscriptomeIndex
from dsforge.controller.design_task_parallel import ParallelDesignTask
from dsforge.controller.design_task import DesignConfig
from dsforge.database.manager import DatabaseManager

def test_parallel_siRNA():
    print("=" * 60)
    print("Test: Parallel siRNA Design (4 cores)")
    print("=" * 60)

    transcriptome = TranscriptomeIndex()
    transcriptome.sequences = {
        "test_gene": "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC",
    }
    transcriptome._compute_stats()

    config = DesignConfig(
        mode="siRNA",
        enabled_rules=["consensus", "reynolds"],
        gc_min=20, gc_max=80,
    )

    task = ParallelDesignTask(
        db_manager=DatabaseManager(":memory:"),
        n_cores=4,
        batch_size=10,
    )

    def progress_cb(step, pct):
        if int(pct) % 25 == 0:
            print(f"  [{pct:5.1f}%] {step}")

    import time
    start = time.time()
    result = task.run_parallel(
        transcriptome=transcriptome,
        target_seq_id="test_gene",
        config=config,
        progress_callback=progress_cb,
    )
    elapsed = time.time() - start

    print(f"\nCompleted in {elapsed:.2f}s")
    print(f"Candidates: {result['summary']['total_candidates']}")
    print(f"Passed: {result['summary']['passed_candidates']}")
    print("✓ Parallel siRNA OK\n")


def test_parallel_long_dsRNA():
    print("=" * 60)
    print("Test: Parallel Long dsRNA Design (4 cores)")
    print("=" * 60)

    long_seq = "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC"

    transcriptome = TranscriptomeIndex()
    transcriptome.sequences = {"test_gene_long": long_seq}
    transcriptome._compute_stats()

    config = DesignConfig(
        mode="long_dsRNA",
        length_min=50, length_max=80,
        enabled_rules=["consensus", "reynolds", "ui_tei"],
        gc_min=20, gc_max=80,
    )

    task = ParallelDesignTask(
        db_manager=DatabaseManager(":memory:"),
        n_cores=4,
        batch_size=50,
    )

    def progress_cb(step, pct):
        if int(pct) % 25 == 0:
            print(f"  [{pct:5.1f}%] {step}")

    import time
    start = time.time()
    result = task.run_parallel(
        transcriptome=transcriptome,
        target_seq_id="test_gene_long",
        config=config,
        progress_callback=progress_cb,
    )
    elapsed = time.time() - start

    print(f"\nCompleted in {elapsed:.2f}s")
    print(f"Regions: {result['summary']['total_candidates']}")
    print(f"Passed: {result['summary']['passed_candidates']}")
    print("✓ Parallel long dsRNA OK\n")


if __name__ == "__main__":
    test_parallel_siRNA()
    test_parallel_long_dsRNA()
    print("=" * 60)
    print("ALL PARALLEL TESTS PASSED ✓")
    print("=" * 60)
