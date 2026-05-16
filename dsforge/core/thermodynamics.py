"""
ViennaRNA 热力学计算封装。

Windows 发布包会随 PyInstaller 一起打入官方 ViennaRNA CLI；RNAup 精筛
只有在实际调用到 `RNAup.exe` 时才标记为 `RNAup-cli`。开发环境缺少 CLI
时仍可返回明确标记的 RNAduplex fallback，便于测试和兼容旧数据。
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional, Dict


class ThermodynamicsCalculator:
    """ViennaRNA 热力学计算器"""

    def __init__(self):
        self._check_vienna()

    def _check_vienna(self):
        """验证 ViennaRNA 是否可用"""
        self.rnaup_executable = self._find_executable("RNAup")
        try:
            import RNA

            self.RNA = RNA
            self.available = True
        except ImportError:
            self.RNA = None
            self.available = False
            print("[Warning] ViennaRNA (RNA module) not available. Thermodynamics disabled.")

    def _find_executable(self, name: str) -> Optional[str]:
        """Find a ViennaRNA CLI in PATH, app folder, PyInstaller temp dir, or installer dirs."""
        candidates = []
        found = shutil.which(name)
        if found:
            return found
        suffixes = [name, f"{name}.exe"]
        for base in [
            Path.cwd(),
            Path(getattr(sys, "_MEIPASS", "")) if getattr(sys, "_MEIPASS", None) else None,
            Path(sys.executable).resolve().parent if sys.executable else None,
            Path(__file__).resolve().parent,
        ]:
            if base:
                candidates.extend(base / suffix for suffix in suffixes)
        vienna_dir = os.environ.get("VIENNA_DLL_DIR")
        if vienna_dir:
            candidates.extend(Path(vienna_dir) / suffix for suffix in suffixes)
        if os.name == "nt":
            for base in [
                Path(r"C:\Program Files\ViennaRNA Package"),
                Path(r"C:\Program Files (x86)\ViennaRNA Package"),
                Path(r"C:\Program Files\ViennaRNA"),
                Path(r"C:\Program Files (x86)\ViennaRNA"),
            ]:
                candidates.extend(base / suffix for suffix in suffixes)
        for candidate in candidates:
            if candidate.exists():
                return str(candidate)
        return None

    def rnaduplex(self, seq1: str, seq2: str) -> Optional[Dict]:
        """
        RNAduplex: 仅计算分子间配对（忽略分子内结构）
        适用于种子区快速预筛

        Returns:
            {
                'dg': float,       # 结合自由能 (kcal/mol)
                'structure': str,  # 配对结构
            }
        """
        if not self.available:
            return None

        try:
            seq1 = seq1.upper().replace("T", "U")
            seq2 = seq2.upper().replace("T", "U")

            # RNAduplex 计算
            duplex = self.RNA.duplexfold(seq1, seq2)
            dg = duplex.energy
            structure = duplex.structure

            return {
                "dg": dg,
                "structure": structure,
            }
        except Exception as e:
            print(f"[RNAduplex Error] {e}")
            return None

    def rnacofold(self, seq1: str, seq2: str) -> Optional[Dict]:
        """
        RNAcofold: 计算两分子杂交结构（含分子内结构）
        适用于精确评估 siRNA-mRNA 配对

        Returns:
            {
                'dg': float,       # 总自由能
                'structure': str,  # 联合结构
            }
        """
        if not self.available:
            return None

        try:
            seq1 = seq1.upper().replace("T", "U")
            seq2 = seq2.upper().replace("T", "U")

            # RNAcofold 计算
            fc = self.RNA.cofold(seq1 + "&" + seq2)
            dg = fc[1]
            structure = fc[0]

            return {
                "dg": dg,
                "structure": structure,
            }
        except Exception as e:
            print(f"[RNAcofold Error] {e}")
            return None

    def rnaup(self, seq1: str, seq2: str) -> Optional[Dict]:
        """
        RNAup: 计算结合自由能，含"打开"mRNA 结构的能量惩罚
        最准确但最慢，建议仅对 top-50 候选运行。

        Note: ViennaRNA 2.7 Python API 中 RNAup 需通过 CLI 调用。
        Windows 打包版要求 `RNAup.exe` 可用；开发环境缺少 CLI 时会明确
        返回 RNAduplex fallback，不把它标记为完整 RNAup 精确结果。

        Returns:
            {
                'dg': float,       # 总结合自由能
                'details': Dict,   # 分解的能量项
            }
        """
        if not self.available:
            return None

        if self.rnaup_executable:
            try:
                seq1 = seq1.upper().replace("T", "U")
                seq2 = seq2.upper().replace("T", "U")
                payload = f">target\n{seq2}\n>query\n{seq1}\n"
                with tempfile.TemporaryDirectory() as tmpdir:
                    proc = subprocess.run(
                        [self.rnaup_executable, "-b"],
                        input=payload,
                        text=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        timeout=30,
                        check=False,
                        cwd=tmpdir,
                    )
                if proc.returncode == 0:
                    match = re.search(r"\((-?\d+(?:\.\d+)?)\s*=", proc.stdout)
                    if not match:
                        match = re.search(r"\((-?\d+(?:\.\d+)?)\)", proc.stdout)
                    if match:
                        return {
                            "dg": float(match.group(1)),
                            "details": {
                                "method": "RNAup-cli",
                                "rnaup_available": True,
                                "stdout": proc.stdout.strip()[:500],
                            },
                        }
            except Exception as e:
                print(f"[RNAup CLI Error] {e}")

        # Fallback: 用 RNAduplex 近似，并明确标记不能当成完整 RNAup。
        fallback = self.rnaduplex(seq1, seq2)
        if fallback:
            return {
                "dg": fallback["dg"],
                "details": {
                    "method": "RNAduplex-fallback-for-RNAup",
                    "rnaup_available": False,
                    "note": "RNAup CLI not available; value is RNAduplex fallback, not full RNAup precision.",
                },
            }
        return None

    def calculate_seed_dg(self, seed: str, target_seq: str) -> Optional[float]:
        """
        计算种子区与靶标序列的 RNAduplex 能量
        用于快速脱靶预筛

        Args:
            seed: 7nt 种子序列
            target_seq: 靶标 mRNA 序列

        Returns:
            最小 ΔG 值（最负 = 结合最强）
        """
        if not self.available or not seed or not target_seq:
            return None

        min_dg = 0.0
        found = False

        # 滑动窗口计算种子区与靶标各位置的杂交能量
        for i in range(len(target_seq) - len(seed) + 1):
            target_sub = target_seq[i : i + len(seed)]
            result = self.rnaduplex(seed, target_sub)
            if result is not None:
                dg = result["dg"]
                if not found or dg < min_dg:
                    min_dg = dg
                    found = True

        return min_dg if found else None
