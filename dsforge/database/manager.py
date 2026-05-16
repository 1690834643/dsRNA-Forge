"""
SQLite 数据库管理器
封装数据库的 CRUD 操作
"""

import sqlite3
import json
from pathlib import Path
from typing import Optional, List, Dict, Any

from dsforge.database.schema import init_database


def default_database_path() -> Path:
    """Return a writable per-user database path for the desktop app."""
    path = Path.home() / ".dsrna_forge"
    path.mkdir(parents=True, exist_ok=True)
    return path / "dsrna_forge.db"


class AutoClosingConnection(sqlite3.Connection):
    """sqlite3 connection that closes after a with-block transaction."""

    def __exit__(self, exc_type, exc_value, traceback):
        result = super().__exit__(exc_type, exc_value, traceback)
        self.close()
        return result


class DatabaseManager:
    """数据库管理器"""

    def __init__(self, db_path: Optional[str] = None):
        self.db_path = Path(db_path) if db_path is not None else default_database_path()
        self._is_memory = str(self.db_path) == ":memory:"
        self._persistent_conn = None
        if self._is_memory:
            self._persistent_conn = sqlite3.connect(":memory:")
            self._persistent_conn.row_factory = sqlite3.Row
            self._persistent_conn.execute("PRAGMA foreign_keys = ON")
            init_database(self._persistent_conn)
        self._init_db()

    def _get_conn(self) -> sqlite3.Connection:
        """获取数据库连接（强制启用外键约束）"""
        if self._is_memory and self._persistent_conn is not None:
            return self._persistent_conn
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(str(self.db_path), factory=AutoClosingConnection)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON")
        # Verify foreign keys are actually enabled (SQLite default is OFF)
        fk_status = conn.execute("PRAGMA foreign_keys").fetchone()[0]
        if not fk_status:
            conn.close()
            raise RuntimeError("Failed to enable SQLite foreign key constraints.")
        return conn

    def _init_db(self):
        """初始化数据库"""
        if not self._is_memory:
            conn = self._get_conn()
            init_database(conn)
            conn.close()

    # === Design Tasks ===

    def create_task(
        self,
        mode: str,
        target_seq_id: Optional[str],
        target_seq: Optional[str],
        params: Dict,
    ) -> int:
        """创建新任务，返回任务 ID"""
        with self._get_conn() as conn:
            cursor = conn.execute(
                """
                INSERT INTO design_tasks (mode, target_seq_id, target_seq, params_json, status)
                VALUES (?, ?, ?, ?, 'pending')
                """,
                (mode, target_seq_id, target_seq, json.dumps(params)),
            )
            conn.commit()
            return cursor.lastrowid

    def update_task_status(self, task_id: int, status: str):
        """更新任务状态"""
        with self._get_conn() as conn:
            if status in ("completed", "cancelled", "failed"):
                conn.execute(
                    """
                    UPDATE design_tasks
                    SET status = ?, completed_at = CURRENT_TIMESTAMP
                    WHERE id = ?
                    """,
                    (status, task_id),
                )
            else:
                conn.execute(
                    "UPDATE design_tasks SET status = ? WHERE id = ?",
                    (status, task_id),
                )
            conn.commit()

    def get_task(self, task_id: int) -> Optional[Dict]:
        """获取任务详情"""
        with self._get_conn() as conn:
            row = conn.execute(
                "SELECT * FROM design_tasks WHERE id = ?", (task_id,)
            ).fetchone()
            if row is None:
                return None
            return dict(row)

    def list_tasks(self, limit: int = 100) -> List[Dict]:
        """列出历史任务"""
        with self._get_conn() as conn:
            rows = conn.execute(
                "SELECT * FROM design_tasks ORDER BY created_at DESC LIMIT ?",
                (limit,),
            ).fetchall()
            return [dict(row) for row in rows]

    def delete_task(self, task_id: int):
        """删除任务（级联删除相关结果）"""
        with self._get_conn() as conn:
            conn.execute("DELETE FROM design_tasks WHERE id = ?", (task_id,))
            conn.commit()

    # === Results ===

    def add_result(self, task_id: int, rank: int, candidate_seq: str,
                   position_start: int, position_end: int,
                   consensus_score: float, passed_filters: int = 0,
                   risk_level: str = "low", risk_score: float = 0,
                   top_risk_targets: str = "", validation_direction: str = "",
                   recommendation_score: float = 0, cluster_id: int = 0,
                   cluster_size: int = 1, alternative_count: int = 0,
                   cluster_span: str = "", explanation_json: str = "{}",
                   validation_hits_json: str = "[]", primers_json: str = "{}",
                   rnaup_json: str = "{}", sgrna_json: str = "{}",
                   region_map: str = "") -> int:
        """添加结果记录"""
        with self._get_conn() as conn:
            cursor = conn.execute(
                """
                INSERT INTO results (
                    task_id, rank, candidate_seq, position_start, position_end,
                    consensus_score, passed_filters, risk_level, risk_score,
                    top_risk_targets, validation_direction, recommendation_score,
                    cluster_id, cluster_size, alternative_count, cluster_span,
                    explanation_json, validation_hits_json, primers_json, rnaup_json,
                    sgrna_json, region_map
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    task_id,
                    rank,
                    candidate_seq,
                    position_start,
                    position_end,
                    consensus_score,
                    passed_filters,
                    risk_level,
                    risk_score,
                    top_risk_targets,
                    validation_direction,
                    recommendation_score,
                    cluster_id,
                    cluster_size,
                    alternative_count,
                    cluster_span,
                    explanation_json,
                    validation_hits_json,
                    primers_json,
                    rnaup_json,
                    sgrna_json,
                    region_map,
                ),
            )
            conn.commit()
            return cursor.lastrowid

    def get_results(self, task_id: int) -> List[Dict]:
        """获取任务的所有结果"""
        with self._get_conn() as conn:
            rows = conn.execute(
                "SELECT * FROM results WHERE task_id = ? ORDER BY rank",
                (task_id,),
            ).fetchall()
            results = []
            for row in rows:
                item = dict(row)
                for column, fallback in [
                    ("explanation_json", {}),
                    ("validation_hits_json", []),
                    ("primers_json", {}),
                    ("rnaup_json", {}),
                    ("sgrna_json", {}),
                ]:
                    try:
                        parsed = json.loads(item.get(column) or json.dumps(fallback))
                    except (json.JSONDecodeError, TypeError):
                        parsed = fallback
                    if column == "explanation_json":
                        item["explanation"] = parsed
                    elif column == "validation_hits_json":
                        item["validation_hits"] = parsed
                    elif column == "primers_json":
                        item["primers"] = parsed
                    elif column == "rnaup_json":
                        item["rnaup"] = parsed
                    elif column == "sgrna_json":
                        item["sgrna"] = parsed
                item["sequence"] = item.get("candidate_seq", "")
                item["position"] = f"{item.get('position_start', '')}-{item.get('position_end', '')}"
                item["passed"] = bool(item.get("passed_filters", False))
                results.append(item)
            return results

    # === Rule Scores ===

    def add_rule_score(self, result_id: int, rule_name: str, score: float,
                       passed: bool, violations: Optional[List[str]] = None):
        """添加规则评分"""
        with self._get_conn() as conn:
            conn.execute(
                """
                INSERT INTO rule_scores (result_id, rule_name, score, passed, violations)
                VALUES (?, ?, ?, ?, ?)
                """,
                (result_id, rule_name, score, int(passed), json.dumps(violations or [])),
            )
            conn.commit()

    def get_rule_scores(self, result_id: int) -> List[Dict]:
        """获取结果的规则评分"""
        with self._get_conn() as conn:
            rows = conn.execute(
                "SELECT * FROM rule_scores WHERE result_id = ?",
                (result_id,),
            ).fetchall()
            return [dict(row) for row in rows]

    # === Thermodynamics ===

    def add_thermodynamics(self, result_id: int, on_target_dg: Optional[float],
                           seed_matches_count: int = 0,
                           high_risk_off_targets: int = 0,
                           rnaup_dg: Optional[float] = None):
        """添加热力学评估"""
        with self._get_conn() as conn:
            conn.execute(
                """
                INSERT INTO thermodynamics (result_id, on_target_dg, seed_matches_count, high_risk_off_targets, rnaup_dg)
                VALUES (?, ?, ?, ?, ?)
                """,
                (result_id, on_target_dg, seed_matches_count, high_risk_off_targets, rnaup_dg),
            )
            conn.commit()

    # === Pool Details ===

    def add_pool_detail(self, result_id: int, dicer_product_seq: str,
                        cut_position: Optional[int], product_score: Optional[float]):
        """添加 Pool 明细"""
        with self._get_conn() as conn:
            conn.execute(
                """
                INSERT INTO pool_details (result_id, dicer_product_seq, cut_position, product_score)
                VALUES (?, ?, ?, ?)
                """,
                (result_id, dicer_product_seq, cut_position, product_score),
            )
            conn.commit()
