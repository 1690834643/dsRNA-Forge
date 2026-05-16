#!/usr/bin/env python3
"""
dsRNA-Forge — 昆虫/植物长 dsRNA 设计工具
入口文件
"""

import sys
import os
import multiprocessing as mp

# 确保项目根目录在路径中
PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

# 跨平台 ViennaRNA 加载（必须在 import RNA 之前）
from dsforge.utils.vienna_loader import setup_vienna_path

setup_vienna_path()


def run_runtime_check() -> int:
    """Verify that the packaged runtime can import GUI/science deps and design."""
    errors = []

    try:
        import PyQt6  # noqa: F401
    except Exception as exc:
        errors.append(f"PyQt6 import failed: {exc}")

    try:
        import RNA
        print(f"ViennaRNA: {getattr(RNA, '__version__', 'available')}")
    except Exception as exc:
        message = f"ViennaRNA import failed: {exc}"
        if os.name == "nt" or getattr(sys, "frozen", False):
            errors.append(message)
        else:
            print(f"WARNING: {message}")

    try:
        from dsforge.core.sequence import TranscriptomeIndex
        from dsforge.controller.design_task import DesignConfig, DesignTask
        from dsforge.database.manager import DatabaseManager
        from dsforge.core.thermodynamics import ThermodynamicsCalculator

        transcriptome = TranscriptomeIndex()
        transcriptome.sequences = {
            "target_gene": "AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGC",
            "background_gene": "UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU",
        }
        transcriptome._compute_stats()
        config = DesignConfig(
            mode="siRNA",
            gc_min=20,
            gc_max=80,
            enabled_rules=["consensus", "reynolds"],
            n_cores=1,
        )
        result = DesignTask(DatabaseManager(":memory:")).run(
            transcriptome=transcriptome,
            target_seq_id="target_gene",
            config=config,
        )
        if result["summary"]["passed_candidates"] <= 0:
            errors.append(f"Runtime design produced no passing candidates: {result['summary']}")
        else:
            print(f"Runtime design OK: {result['summary']}")

        thermo = ThermodynamicsCalculator()
        rnaup = thermo.rnaup("AUGCAUGCAUGCAUGCAUGC", "GCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAU")
        rnaup_method = ((rnaup or {}).get("details") or {}).get("method")
        print(f"RNAup method: {rnaup_method or 'unavailable'}")
        if (os.name == "nt" or getattr(sys, "frozen", False)) and rnaup_method != "RNAup-cli":
            errors.append(f"RNAup CLI is required in Windows packaged runtime, got: {rnaup_method}")

        sgrna_ref = TranscriptomeIndex()
        sgrna_ref.sequences = {
            "target_gene": "A" * 25 + "G" * 20 + "AGG" + "T" * 40,
            "background_gene": "C" * 100,
        }
        sgrna_ref._compute_stats()
        sgrna_result = DesignTask(DatabaseManager(":memory:")).run(
            transcriptome=sgrna_ref,
            target_seq_id="target_gene",
            config=DesignConfig(mode="sgRNA", n_cores=1),
        )
        if sgrna_result["summary"]["total_candidates"] <= 0:
            errors.append(f"Runtime sgRNA design produced no candidates: {sgrna_result['summary']}")
        else:
            print(f"Runtime sgRNA OK: {sgrna_result['summary']}")
    except Exception as exc:
        errors.append(f"Runtime design failed: {exc}")

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print("Runtime check passed")
    return 0


def main():
    """主入口"""
    import argparse

    parser = argparse.ArgumentParser(
        description="dsRNA-Forge: dsRNA/siRNA Design Tool for Insects and Plants"
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        default=True,
        help="Launch GUI (default)",
    )
    parser.add_argument(
        "--cli",
        action="store_true",
        help="Run in CLI mode",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="dsRNA-Forge 0.1.0",
    )
    parser.add_argument(
        "--check-runtime",
        action="store_true",
        help="Verify packaged dependencies and run a tiny design, then exit.",
    )
    args = parser.parse_args()

    if args.check_runtime:
        sys.exit(run_runtime_check())

    if args.cli:
        print("CLI mode not yet implemented. Use --gui or run without arguments.")
        sys.exit(1)

    # GUI 模式
    from PyQt6.QtWidgets import QApplication
    from dsforge.gui.main_window import MainWindow

    app = QApplication(sys.argv)
    app.setApplicationName("dsRNA-Forge")
    app.setApplicationVersion("0.1.0")

    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    mp.freeze_support()
    main()
