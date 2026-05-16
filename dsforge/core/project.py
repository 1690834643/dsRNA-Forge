"""
Project save/load helpers for dsRNA-Forge.
"""

import json
import time
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Dict


PROJECT_VERSION = 1


def _jsonable(value: Any):
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    return value


def save_project_file(filepath: str, payload: Dict):
    """Save a portable project JSON file."""
    path = Path(filepath)
    data = {
        "version": PROJECT_VERSION,
        "saved_at": time.time(),
        **_jsonable(payload),
    }
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def load_project_file(filepath: str) -> Dict:
    """Load and validate a project JSON file."""
    path = Path(filepath)
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("version") != PROJECT_VERSION:
        raise ValueError(f"Unsupported dsRNA-Forge project version: {data.get('version')}")
    return data
