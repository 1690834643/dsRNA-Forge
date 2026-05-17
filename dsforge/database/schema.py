"""
SQLite 数据库 Schema
规范化设计，支持三种设计模式的结果存储
"""

import re

SCHEMA_SQL = """
-- 设计任务表
CREATE TABLE IF NOT EXISTS design_tasks (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    mode TEXT NOT NULL CHECK(mode IN ('siRNA', 'DsiRNA', 'long_dsRNA', 'sgRNA')),
    target_seq_id TEXT,
    target_seq TEXT,
    params_json TEXT NOT NULL,
    status TEXT NOT NULL DEFAULT 'pending' CHECK(status IN ('pending', 'running', 'completed', 'cancelled', 'failed')),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    completed_at TIMESTAMP
);

-- 通用结果表（三种模式共用）
CREATE TABLE IF NOT EXISTS results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    task_id INTEGER NOT NULL,
    rank INTEGER NOT NULL,
    candidate_seq TEXT NOT NULL,
    position_start INTEGER,
    position_end INTEGER,
    consensus_score REAL,
    passed_filters INTEGER DEFAULT 0,
    FOREIGN KEY (task_id) REFERENCES design_tasks(id) ON DELETE CASCADE
);

-- 规则评分明细（每条 result 对应多条）
CREATE TABLE IF NOT EXISTS rule_scores (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    result_id INTEGER NOT NULL,
    rule_name TEXT NOT NULL,
    score REAL,
    passed INTEGER DEFAULT 0,
    violations TEXT,
    FOREIGN KEY (result_id) REFERENCES results(id) ON DELETE CASCADE
);

-- 热力学评估（脱靶摘要）
CREATE TABLE IF NOT EXISTS thermodynamics (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    result_id INTEGER NOT NULL,
    on_target_dg REAL,
    seed_matches_count INTEGER DEFAULT 0,
    high_risk_off_targets INTEGER DEFAULT 0,
    rnaup_dg REAL,
    FOREIGN KEY (result_id) REFERENCES results(id) ON DELETE CASCADE
);

-- 长 dsRNA 特有：Pool 明细
CREATE TABLE IF NOT EXISTS pool_details (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    result_id INTEGER NOT NULL,
    dicer_product_seq TEXT NOT NULL,
    cut_position INTEGER,
    product_score REAL,
    FOREIGN KEY (result_id) REFERENCES results(id) ON DELETE CASCADE
);

-- 索引优化
CREATE INDEX IF NOT EXISTS idx_results_task_id ON results(task_id);
CREATE INDEX IF NOT EXISTS idx_rule_scores_result_id ON rule_scores(result_id);
CREATE INDEX IF NOT EXISTS idx_thermodynamics_result_id ON thermodynamics(result_id);
CREATE INDEX IF NOT EXISTS idx_pool_details_result_id ON pool_details(result_id);
"""

RESULTS_TABLE_SQL = """
CREATE TABLE results (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    task_id INTEGER NOT NULL,
    rank INTEGER NOT NULL,
    candidate_seq TEXT NOT NULL,
    position_start INTEGER,
    position_end INTEGER,
    consensus_score REAL,
    passed_filters INTEGER DEFAULT 0,
    risk_level TEXT DEFAULT 'low',
    risk_score REAL DEFAULT 0,
    top_risk_targets TEXT DEFAULT '',
    validation_direction TEXT DEFAULT '',
    recommendation_score REAL DEFAULT 0,
    cluster_id INTEGER DEFAULT 0,
    cluster_size INTEGER DEFAULT 1,
    alternative_count INTEGER DEFAULT 0,
    cluster_span TEXT DEFAULT '',
    off_target_json TEXT DEFAULT '{}',
    explanation_json TEXT DEFAULT '{}',
    validation_hits_json TEXT DEFAULT '[]',
    primers_json TEXT DEFAULT '{}',
    rnaup_json TEXT DEFAULT '{}',
    sgrna_json TEXT DEFAULT '{}',
    region_map TEXT DEFAULT '',
    FOREIGN KEY (task_id) REFERENCES design_tasks(id) ON DELETE CASCADE
);
"""

RESULTS_COLUMNS = [
    "id",
    "task_id",
    "rank",
    "candidate_seq",
    "position_start",
    "position_end",
    "consensus_score",
    "passed_filters",
    "risk_level",
    "risk_score",
    "top_risk_targets",
    "validation_direction",
    "recommendation_score",
    "cluster_id",
    "cluster_size",
    "alternative_count",
    "cluster_span",
    "off_target_json",
    "explanation_json",
    "validation_hits_json",
    "primers_json",
    "rnaup_json",
    "sgrna_json",
    "region_map",
]


_TABLE_COL_DEF_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _ensure_column(conn, table: str, column: str, definition: str):
    # Validate identifiers to prevent SQL injection
    for identifier in (table, column):
        if not _TABLE_COL_DEF_RE.match(identifier):
            raise ValueError(f"Invalid SQL identifier: {identifier!r}")
    columns = [row[1] for row in conn.execute(f"PRAGMA table_info({table})").fetchall()]
    if column not in columns:
        conn.execute(f"ALTER TABLE {table} ADD COLUMN {column} {definition}")


def _set_foreign_keys(conn, enabled: bool):
    conn.commit()
    conn.execute(f"PRAGMA foreign_keys = {'ON' if enabled else 'OFF'}")


def _recreate_design_tasks_with_sgrna(conn):
    previous_fk = bool(conn.execute("PRAGMA foreign_keys").fetchone()[0])
    _set_foreign_keys(conn, False)
    conn.execute("PRAGMA legacy_alter_table = ON")
    conn.executescript(
        """
        ALTER TABLE design_tasks RENAME TO design_tasks_old;
        CREATE TABLE design_tasks (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            mode TEXT NOT NULL CHECK(mode IN ('siRNA', 'DsiRNA', 'long_dsRNA', 'sgRNA')),
            target_seq_id TEXT,
            target_seq TEXT,
            params_json TEXT NOT NULL,
            status TEXT NOT NULL DEFAULT 'pending' CHECK(status IN ('pending', 'running', 'completed', 'cancelled', 'failed')),
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            completed_at TIMESTAMP
        );
        INSERT INTO design_tasks (id, mode, target_seq_id, target_seq, params_json, status, created_at, completed_at)
        SELECT id, mode, target_seq_id, target_seq, params_json, status, created_at, completed_at
        FROM design_tasks_old;
        DROP TABLE design_tasks_old;
        """
    )
    conn.execute("PRAGMA legacy_alter_table = OFF")
    conn.commit()
    if previous_fk:
        _set_foreign_keys(conn, True)


def _rebuild_results_table(conn):
    """Repair results foreign keys after legacy design_tasks migrations."""
    existing_columns = [
        row[1] for row in conn.execute("PRAGMA table_info(results)").fetchall()
    ]
    if not existing_columns:
        return

    copy_columns = [column for column in RESULTS_COLUMNS if column in existing_columns]
    previous_fk = bool(conn.execute("PRAGMA foreign_keys").fetchone()[0])
    _set_foreign_keys(conn, False)
    conn.execute("PRAGMA legacy_alter_table = ON")
    conn.execute("ALTER TABLE results RENAME TO results_old")
    conn.execute(RESULTS_TABLE_SQL)
    if copy_columns:
        column_sql = ", ".join(copy_columns)
        conn.execute(
            f"INSERT INTO results ({column_sql}) SELECT {column_sql} FROM results_old"
        )
    conn.execute("DROP TABLE results_old")
    conn.execute("PRAGMA legacy_alter_table = OFF")
    conn.commit()
    if previous_fk:
        _set_foreign_keys(conn, True)


def _results_fk_targets(conn):
    return [row[2] for row in conn.execute("PRAGMA foreign_key_list(results)").fetchall()]


def init_database(conn):
    """初始化数据库，创建所有表"""
    conn.executescript(SCHEMA_SQL)
    task_sql = conn.execute(
        "SELECT sql FROM sqlite_master WHERE type='table' AND name='design_tasks'"
    ).fetchone()
    if task_sql and "sgRNA" not in (task_sql[0] or ""):
        _recreate_design_tasks_with_sgrna(conn)
    _ensure_column(conn, "results", "risk_level", "TEXT DEFAULT 'low'")
    _ensure_column(conn, "results", "risk_score", "REAL DEFAULT 0")
    _ensure_column(conn, "results", "top_risk_targets", "TEXT DEFAULT ''")
    _ensure_column(conn, "results", "validation_direction", "TEXT DEFAULT ''")
    _ensure_column(conn, "results", "recommendation_score", "REAL DEFAULT 0")
    _ensure_column(conn, "results", "cluster_id", "INTEGER DEFAULT 0")
    _ensure_column(conn, "results", "cluster_size", "INTEGER DEFAULT 1")
    _ensure_column(conn, "results", "alternative_count", "INTEGER DEFAULT 0")
    _ensure_column(conn, "results", "cluster_span", "TEXT DEFAULT ''")
    _ensure_column(conn, "results", "off_target_json", "TEXT DEFAULT '{}'")
    _ensure_column(conn, "results", "explanation_json", "TEXT DEFAULT '{}'")
    _ensure_column(conn, "results", "validation_hits_json", "TEXT DEFAULT '[]'")
    _ensure_column(conn, "results", "primers_json", "TEXT DEFAULT '{}'")
    _ensure_column(conn, "results", "rnaup_json", "TEXT DEFAULT '{}'")
    _ensure_column(conn, "results", "sgrna_json", "TEXT DEFAULT '{}'")
    _ensure_column(conn, "results", "region_map", "TEXT DEFAULT ''")
    if _results_fk_targets(conn) != ["design_tasks"]:
        _rebuild_results_table(conn)
    conn.commit()
