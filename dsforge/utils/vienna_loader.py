"""
跨平台 ViennaRNA 加载封装
处理 Linux / Windows 的 RNA 模块导入差异
"""

import os
import sys
import platform


def setup_vienna_path():
    """
    配置 ViennaRNA 库路径。
    必须在 import RNA 之前调用。
    """
    system = platform.system()

    if system == "Windows":
        _setup_windows()
    else:
        _setup_linux()


def _setup_windows():
    """Windows: 通过 os.add_dll_directory 添加 RNA.dll 路径"""
    # 尝试多个可能的路径
    possible_paths = [
        # PyInstaller one-file extraction directory
        getattr(sys, "_MEIPASS", ""),
        # PyInstaller 打包后的同级目录
        os.path.dirname(sys.executable),
        # 开发环境下的 ViennaRNA 安装路径
        r"C:\Program Files\ViennaRNA Package",
        r"C:\Program Files (x86)\ViennaRNA Package",
        r"C:\Program Files\ViennaRNA",
        r"C:\Program Files (x86)\ViennaRNA",
        # 当前目录
        os.getcwd(),
    ]

    # 也检查环境变量
    if "VIENNARNA_PATH" in os.environ:
        possible_paths.insert(0, os.environ["VIENNARNA_PATH"])
    if "VIENNA_DLL_DIR" in os.environ:
        possible_paths.insert(0, os.environ["VIENNA_DLL_DIR"])

    for path in possible_paths:
        dll_path = os.path.join(path, "RNA.dll")
        if os.path.exists(dll_path):
            try:
                os.add_dll_directory(path)
                print(f"[vienna_loader] Added DLL directory: {path}")
                return
            except Exception as e:
                print(f"[vienna_loader] Failed to add {path}: {e}")

    # The pip wheel can load its bundled extension/DLLs without RNA.dll being
    # visible in standard install directories. Do not warn if import already works.
    try:
        import RNA  # noqa: F401
        return
    except Exception:
        print("[vienna_loader] Warning: RNA.dll not found in standard paths.")
        print("[vienna_loader] Set VIENNARNA_PATH environment variable if needed.")


def _setup_linux():
    """Linux: 通常 pip 安装的 ViennaRNA 可直接 import"""
    # pip 安装的 ViennaRNA 通常不需要额外路径设置
    # 如果需要系统级 libRNA，可以在这里设置 LD_LIBRARY_PATH
    pass


def get_vienna_version() -> str:
    """获取 ViennaRNA 版本号，用于验证安装"""
    try:
        import RNA

        return RNA.__version__
    except ImportError:
        return "NOT_INSTALLED"


def check_vienna_available() -> bool:
    """检查 ViennaRNA 是否可用"""
    try:
        import RNA

        # 做一个简单的测试调用
        seq = "GGGAAACCC"
        fc = RNA.fold_compound(seq)
        return True
    except Exception:
        return False


if __name__ == "__main__":
    setup_vienna_path()
    version = get_vienna_version()
    available = check_vienna_available()
    print(f"ViennaRNA version: {version}")
    print(f"ViennaRNA available: {available}")
