"""
序列处理模块
- FASTA 加载与索引
- 序列工具函数
- 候选生成器（滑动窗口）
"""

import os
import hashlib
import json
import re
import time
from typing import Any, Dict, List, Tuple, Iterator, Optional
from pathlib import Path


# 默认缓存目录
DEFAULT_CACHE_DIR = Path.home() / ".cache" / "dsrna_forge"
VALID_RNA_BASES = set("AUGC")
AMBIGUOUS_RNA_BASES = set("NRYMKSWHBVD")


def _resolve_cache_dir(cache_dir: Optional[Path] = None) -> Path:
    return Path(cache_dir) if cache_dir is not None else DEFAULT_CACHE_DIR


def _manifest_path(cache_dir: Optional[Path] = None) -> Path:
    return _resolve_cache_dir(cache_dir) / "transcriptomes.json"


def _read_manifest(cache_dir: Optional[Path] = None) -> List[Dict]:
    path = _manifest_path(cache_dir)
    if not path.exists():
        return []
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, list) else []
    except Exception:
        return []


def _write_manifest(entries: List[Dict], cache_dir: Optional[Path] = None):
    cache_path = _resolve_cache_dir(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    with open(_manifest_path(cache_path), "w", encoding="utf-8") as f:
        json.dump(entries, f, ensure_ascii=False, indent=2)


def normalize_sequence(sequence: str, strict: bool = False) -> str:
    """Normalize pasted or FASTA sequence text to uppercase RNA letters.

    Args:
        sequence: Input sequence string.
        strict: If True, raise ValueError when non-AUGCT characters are present.
    """
    cleaned = re.sub(r"\s+", "", sequence).upper().replace("T", "U")
    if strict:
        invalid = set(cleaned) - VALID_RNA_BASES
        if invalid:
            raise ValueError(
                f"Sequence contains invalid characters: {sorted(invalid)}. "
                f"Only A, U, G, C (and T) are allowed."
            )
    return cleaned


def make_safe_sequence_id(raw_id: str, fallback: str = "custom_target") -> str:
    """Convert a FASTA header into a stable simple sequence ID."""
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", raw_id.strip()).strip("_")
    return cleaned or fallback


def parse_fasta_records(filepath: str) -> List[Tuple[str, str]]:
    """Parse FASTA records and return normalized RNA sequences."""
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {filepath}")

    records: List[Tuple[str, str]] = []
    current_id = None
    current_seq = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((make_safe_sequence_id(current_id), normalize_sequence("".join(current_seq))))
                current_id = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)

    if current_id is not None:
        records.append((make_safe_sequence_id(current_id), normalize_sequence("".join(current_seq))))

    if not records:
        raise ValueError("No FASTA records found")
    return records


def load_first_fasta_record(filepath: str, fallback_id: str = "custom_target") -> Tuple[str, str]:
    """Load the first FASTA record for custom target design."""
    record_id, sequence = parse_fasta_records(filepath)[0]
    return record_id or fallback_id, sequence


def clone_with_custom_sequence(
    index: "TranscriptomeIndex",
    seq_id: str,
    sequence: str,
    fallback_id: str = "custom_target",
) -> "TranscriptomeIndex":
    """Return a copy of a transcriptome with one custom target added."""
    cloned = TranscriptomeIndex(cache_dir=getattr(index, "cache_dir", None))
    cloned.sequences = dict(index.sequences)
    safe_id = make_safe_sequence_id(seq_id, fallback_id)
    if safe_id in cloned.sequences:
        safe_id = make_safe_sequence_id(f"{safe_id}_custom", fallback_id)
    cloned.sequences[safe_id] = normalize_sequence(sequence)
    cloned.source_file = index.source_file
    digest = hashlib.sha256()
    digest.update(str(index.source_hash or "no_source_hash").encode())
    digest.update(safe_id.encode())
    digest.update(cloned.sequences[safe_id].encode())
    cloned.source_hash = f"{index.source_hash or 'custom'}_target_{digest.hexdigest()[:16]}"
    cloned.cache_path = index.cache_path
    cloned._compute_stats()
    return cloned


def merge_background_transcriptomes(
    primary: "TranscriptomeIndex",
    backgrounds: List[Tuple[str, "TranscriptomeIndex"]],
) -> "TranscriptomeIndex":
    """Return a transcriptome copy with additional off-target backgrounds.

    Background sequence IDs are prefixed with a safe source label so risk reports
    still show where a potential off-target came from.
    """
    merged = TranscriptomeIndex(cache_dir=getattr(primary, "cache_dir", None))
    merged.sequences = dict(primary.sequences)
    merged.source_file = primary.source_file
    merged.cache_path = primary.cache_path

    digest = hashlib.sha256()
    digest.update(str(primary.source_hash or "no_source_hash").encode())
    for seq_id, sequence in sorted(primary.sequences.items()):
        digest.update(seq_id.encode())
        digest.update(str(len(sequence)).encode())

    for label, background in backgrounds or []:
        prefix = make_safe_sequence_id(label, "background")
        digest.update(prefix.encode())
        digest.update(str(getattr(background, "source_hash", "") or "").encode())
        for seq_id, sequence in background.sequences.items():
            normalized = normalize_sequence(sequence)
            digest.update(seq_id.encode())
            digest.update(normalized.encode())
            merged_id = f"{prefix}|{seq_id}"
            suffix = 2
            while merged_id in merged.sequences:
                merged_id = f"{prefix}|{seq_id}_{suffix}"
                suffix += 1
            merged.sequences[merged_id] = normalized

    if backgrounds:
        merged.source_hash = f"{primary.source_hash or 'merged'}_bg_{digest.hexdigest()[:16]}"
    else:
        merged.source_hash = primary.source_hash
    merged._compute_stats()
    return merged


class TranscriptomeIndex:
    """
    转录组索引
    将 FASTA 文件加载为内存中的字典，支持快速查询
    支持索引缓存，加速重复加载
    """

    def __init__(self, cache_dir: Optional[Path] = None):
        self.sequences: Dict[str, str] = {}  # seq_id -> sequence
        self.stats: Dict[str, Any] = {}
        self.source_file: Optional[str] = None
        self.source_hash: Optional[str] = None
        self.cache_dir = _resolve_cache_dir(cache_dir)
        self.cache_path: Optional[str] = None

    def _compute_file_hash(self, filepath: Path) -> str:
        """计算文件完整 SHA256 哈希用于缓存验证"""
        h = hashlib.sha256()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                h.update(chunk)
        return h.hexdigest()[:16]

    def _get_cache_path(self, filepath: Path) -> Path:
        """获取缓存文件路径"""
        cache_dir = self.cache_dir
        cache_dir.mkdir(parents=True, exist_ok=True)
        file_hash = self._compute_file_hash(filepath)
        return cache_dir / f"{filepath.stem}_{file_hash}.idx"

    def load_fasta(self, filepath: str, use_cache: bool = True) -> "TranscriptomeIndex":
        """
        加载 FASTA 文件

        Args:
            filepath: FASTA 文件路径
            use_cache: 是否使用缓存加速（默认开启）
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"FASTA file not found: {filepath}")

        self.source_file = str(filepath)
        self.source_hash = self._compute_file_hash(filepath)
        cache_path = self._get_cache_path(filepath)
        self.cache_path = str(cache_path)

        # 尝试从缓存加载
        if use_cache:
            if cache_path.exists():
                try:
                    self._load_from_cache(cache_path)
                    self.source_file = str(filepath)
                    self.source_hash = self.source_hash or self._compute_file_hash(filepath)
                    self.cache_path = str(cache_path)
                    self._register_saved(filepath, cache_path)
                    return self
                except Exception:
                    pass  # 缓存损坏，重新生成

        # 从 FASTA 解析
        self.sequences = {}
        current_id = None
        current_seq = []

        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_id is not None:
                        self.sequences[current_id] = normalize_sequence("".join(current_seq))
                    current_id = line[1:].split()[0]  # 取第一个空格前的 ID
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_id is not None:
                self.sequences[current_id] = normalize_sequence("".join(current_seq))

        if not self.sequences:
            raise ValueError("No FASTA records found")

        self._compute_stats()

        # 保存缓存
        if use_cache:
            self._save_to_cache(cache_path)
            self._register_saved(filepath, cache_path)

        return self

    def _load_from_cache(self, cache_path: Path):
        """从缓存加载"""
        with open(cache_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if data.get("schema_version") != 2:
            raise ValueError("Unsupported transcriptome cache schema")
        self.sequences = data["sequences"]
        self.stats = data["stats"]
        self.source_file = data.get("source_file", self.source_file)
        self.source_hash = data.get("source_hash", self.source_hash)
        self.cache_path = str(cache_path)

    def _save_to_cache(self, cache_path: Path):
        """保存到缓存"""
        data = {
            "sequences": self.sequences,
            "stats": self.stats,
            "source_file": self.source_file,
            "source_hash": self.source_hash,
            "schema_version": 2,
        }
        with open(cache_path, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False)

    def _register_saved(self, filepath: Path, cache_path: Path):
        """Record this transcriptome in the saved transcriptome manifest."""
        if not self.source_hash:
            return
        entries = [e for e in _read_manifest(self.cache_dir) if e.get("key") != self.source_hash]
        entries.insert(0, {
            "key": self.source_hash,
            "name": filepath.stem,
            "source_file": str(filepath),
            "cache_path": str(cache_path),
            "stats": self.stats,
            "saved_at": time.time(),
        })
        _write_manifest(entries[:50], self.cache_dir)

    @classmethod
    def list_saved(cls, cache_dir: Optional[Path] = None) -> List[Dict]:
        """List transcriptomes that can be loaded from the local cache."""
        entries = []
        for entry in _read_manifest(cache_dir):
            cache_path = Path(entry.get("cache_path", ""))
            if cache_path.exists():
                entries.append(entry)
        return entries

    @classmethod
    def load_saved(cls, key: str, cache_dir: Optional[Path] = None) -> "TranscriptomeIndex":
        """Load a previously registered transcriptome from its cached index."""
        for entry in cls.list_saved(cache_dir):
            if entry.get("key") == key:
                index = cls(cache_dir=cache_dir)
                index.source_file = entry.get("source_file")
                index.source_hash = entry.get("key")
                index.cache_path = entry.get("cache_path")
                index._load_from_cache(Path(entry["cache_path"]))
                index.source_file = entry.get("source_file")
                index.source_hash = entry.get("key")
                return index
        raise KeyError(f"Saved transcriptome not found: {key}")

    @classmethod
    def rename_saved(cls, key: str, new_name: str, cache_dir: Optional[Path] = None):
        """Rename a saved transcriptome in the local manifest."""
        safe_name = new_name.strip()
        if not safe_name:
            raise ValueError("Saved transcriptome name cannot be empty")
        entries = _read_manifest(cache_dir)
        found = False
        for entry in entries:
            if entry.get("key") == key:
                entry["name"] = safe_name
                found = True
                break
        if not found:
            raise KeyError(f"Saved transcriptome not found: {key}")
        _write_manifest(entries, cache_dir)

    @classmethod
    def delete_saved(
        cls,
        key: str,
        cache_dir: Optional[Path] = None,
        delete_cache: bool = True,
    ):
        """Remove a saved transcriptome entry and optionally its cached index."""
        entries = _read_manifest(cache_dir)
        removed = []
        kept = []
        for entry in entries:
            if entry.get("key") == key:
                removed.append(entry)
            else:
                kept.append(entry)
        if not removed:
            raise KeyError(f"Saved transcriptome not found: {key}")

        if delete_cache:
            for entry in removed:
                cache_path = Path(entry.get("cache_path", ""))
                if cache_path.exists():
                    try:
                        cache_path.unlink()
                    except OSError:
                        pass
                for risk_path in [
                    _resolve_cache_dir(cache_dir) / f"offtarget_{entry.get('key')}_v3_7_16_20_27.json",
                    # Legacy formats (safe to clean up if present)
                    _resolve_cache_dir(cache_dir) / f"offtarget_{entry.get('key')}_7_16_20.pkl",
                    _resolve_cache_dir(cache_dir) / f"offtarget_{entry.get('key')}_v2_7_16_20_27.pkl",
                ]:
                    if risk_path.exists():
                        try:
                            risk_path.unlink()
                        except OSError:
                            pass

        _write_manifest(kept, cache_dir)

    @classmethod
    def clear_cache(cls):
        """清除所有缓存"""
        cache_dir = DEFAULT_CACHE_DIR
        if cache_dir.exists():
            for f in cache_dir.glob("*.idx"):
                f.unlink()
            for f in cache_dir.glob("offtarget_*.json"):
                f.unlink()
            for f in cache_dir.glob("offtarget_*.pkl"):
                f.unlink()  # Legacy cleanup
            manifest = _manifest_path(cache_dir)
            if manifest.exists():
                manifest.unlink()

    @classmethod
    def get_cache_info(cls) -> List[Dict]:
        """获取缓存信息"""
        cache_dir = DEFAULT_CACHE_DIR
        if not cache_dir.exists():
            return []
        result = []
        for f in cache_dir.glob("*.idx"):
            stat = f.stat()
            result.append({
                "file": f.name,
                "size_mb": round(stat.st_size / (1024 * 1024), 2),
                "created": stat.st_mtime,
            })
        return result

    def _compute_stats(self):
        """计算转录组统计信息"""
        if not self.sequences:
            self.stats = {}
            return

        lengths = [len(seq) for seq in self.sequences.values()]
        total_nt = sum(lengths)
        total_gc = sum(seq.count("G") + seq.count("C") for seq in self.sequences.values())

        self.stats = {
            "num_sequences": len(self.sequences),
            "total_nt": total_nt,
            "avg_length": total_nt / len(self.sequences),
            "min_length": min(lengths),
            "max_length": max(lengths),
            "gc_content": (total_gc / total_nt * 100) if total_nt > 0 else 0,
        }

    def get_stats(self) -> Dict:
        """获取统计信息"""
        return self.stats

    def get_sequence(self, seq_id: str) -> Optional[str]:
        """通过 ID 获取序列"""
        return self.sequences.get(seq_id)

    def list_ids(self) -> List[str]:
        """列出所有序列 ID"""
        return list(self.sequences.keys())


def gc_content(sequence: str) -> float:
    """计算 GC 含量百分比"""
    sequence = sequence.upper()
    if not sequence:
        return 0.0
    gc = sequence.count("G") + sequence.count("C")
    return (gc / len(sequence)) * 100


def has_poly_repeat(sequence: str, n: int = 4) -> Tuple[bool, str]:
    """
    检查是否有连续重复的碱基

    Returns:
        (has_repeat, base) — has_repeat 为 True 时 base 是重复的碱基
    """
    sequence = sequence.upper()
    for base in "AUGC":
        if base * n in sequence:
            return True, base
    return False, ""


def generate_candidates(
    target_seq: str,
    mode: str,
    min_len: int,
    max_len: int,
    gc_min: float = 0,
    gc_max: float = 100,
    exclude_poly: int = 4,
    max_candidates: int = 0,
) -> Iterator[Dict]:
    """
    滑动窗口候选生成器

    Args:
        target_seq: 目标序列（RNA，U 代替 T）
        mode: 'siRNA', 'DsiRNA', 'long_dsRNA'
        min_len: 最小长度
        max_len: 最大长度
        gc_min: 最小 GC 含量
        gc_max: 最大 GC 含量
        exclude_poly: 排除连续重复碱基的长度阈值
        max_candidates: 最大候选数上限（0 表示无限制）

    Yields:
        {
            'sequence': str,
            'start': int,
            'end': int,
            'length': int,
            'gc': float,
        }
    """
    target_seq = target_seq.upper().replace("T", "U")
    count = 0

    for length in range(min_len, max_len + 1):
        for start in range(0, len(target_seq) - length + 1):
            if max_candidates > 0 and count >= max_candidates:
                return

            seq = target_seq[start : start + length]
            # Skip windows containing ambiguous or otherwise non-canonical bases.
            if set(seq) - VALID_RNA_BASES:
                continue

            # GC 过滤
            gc = gc_content(seq)
            if not (gc_min <= gc <= gc_max):
                continue

            # 连续重复过滤
            has_repeat, _ = has_poly_repeat(seq, exclude_poly)
            if has_repeat:
                continue

            count += 1
            yield {
                "sequence": seq,
                "start": start,
                "end": start + length,
                "length": length,
                "gc": gc,
            }
