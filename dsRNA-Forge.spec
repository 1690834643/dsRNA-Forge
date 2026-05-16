# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller spec file for dsRNA-Forge
Windows 打包配置

Usage:
    pyinstaller dsRNA-Forge.spec

Requirements:
    - ViennaRNA Windows 版已安装
    - RNA.dll, libgomp-1.dll, vcruntime140.dll 可从 ViennaRNA 安装目录复制
"""

import sys
import os

binaries = []
# Optional external ViennaRNA CLI/DLL directory. The Python wheel provides the
# RNA module; the full Windows installer provides RNAup.exe and CLI DLLs.
vienna_candidates = []
VIENNA_DLL_DIR = os.environ.get("VIENNA_DLL_DIR")
if VIENNA_DLL_DIR:
    vienna_candidates.append(VIENNA_DLL_DIR)
if sys.platform == "win32":
    vienna_candidates.extend([
        r"C:\Program Files\ViennaRNA Package",
        r"C:\Program Files (x86)\ViennaRNA Package",
    ])

vienna_cli_dir = next((path for path in vienna_candidates if path and os.path.exists(os.path.join(path, "RNAup.exe"))), None)
if sys.platform == "win32" and vienna_cli_dir:
    for file_name in os.listdir(vienna_cli_dir):
        lower = file_name.lower()
        if lower.endswith(".dll") or lower in {"rnaup.exe", "rnaduplex.exe", "rnacofold.exe"}:
            binaries.append((os.path.join(vienna_cli_dir, file_name), "."))

a = Analysis(
    ['main.py'],
    pathex=[],
    binaries=binaries,
    datas=[
        ('config.json', '.'),
        ('demo_data', 'demo_data'),
    ],
    hiddenimports=[
        'RNA',
        'dsforge.core.scoring.reynolds',
        'dsforge.core.scoring.consensus',
        'dsforge.core.scoring.ui_tei',
        'dsforge.core.scoring.amarzguioui',
        'dsforge.core.scoring.hsieh',
        'dsforge.core.scoring.jagla',
        'openpyxl',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=None,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=None)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='dsRNA-Forge',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,  # Windows GUI 应用，无控制台窗口
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,  # 可添加 .ico 文件
)
