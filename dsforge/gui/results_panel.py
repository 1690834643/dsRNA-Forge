"""
结果表格面板
- 排序、筛选
- 查看明细
- 导出按钮
"""

from PyQt6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QPushButton,
    QHeaderView,
    QGroupBox,
    QAbstractItemView,
    QTextEdit,
)
from PyQt6.QtCore import Qt


class ResultsPanel(QGroupBox):
    """结果表格面板"""

    def __init__(self):
        super().__init__("Results")
        self.results = []
        self._setup_ui()

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        # 结果表格
        self.table = QTableWidget()
        self.table.setColumnCount(10)
        self.table.setHorizontalHeaderLabels(
            [
                "Rank",
                "Sequence",
                "Position",
                "Consensus Score",
                "Cluster Size",
                "Risk",
                "Risk Score",
                "Top Risk Targets",
                "Validation Direction",
                "Pass",
            ]
        )
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSortingEnabled(True)
        self.table.itemSelectionChanged.connect(self._show_selected_detail)
        layout.addWidget(self.table)

        self.detail_edit = QTextEdit()
        self.detail_edit.setReadOnly(True)
        self.detail_edit.setMaximumHeight(190)
        self.detail_edit.setPlaceholderText(
            "Select a recommendation to inspect ranking reasons, off-target validation, region map and primers/oligos."
        )
        layout.addWidget(self.detail_edit)

        # 操作按钮
        btn_layout = QHBoxLayout()
        self.export_csv_btn = QPushButton("Export CSV")
        self.export_excel_btn = QPushButton("Export Excel")
        self.export_fasta_btn = QPushButton("Export FASTA")
        self.export_report_btn = QPushButton("Export Report")
        self.export_primers_btn = QPushButton("Export Primers")
        self.save_project_btn = QPushButton("Save Project")
        self.open_project_btn = QPushButton("Open Project")

        btn_layout.addWidget(self.export_csv_btn)
        btn_layout.addWidget(self.export_excel_btn)
        btn_layout.addWidget(self.export_fasta_btn)
        btn_layout.addWidget(self.export_report_btn)
        btn_layout.addWidget(self.export_primers_btn)
        btn_layout.addWidget(self.save_project_btn)
        btn_layout.addWidget(self.open_project_btn)
        btn_layout.addStretch()

        layout.addLayout(btn_layout)

    def clear_results(self):
        """清空结果"""
        self.results = []
        self.table.setRowCount(0)
        self.detail_edit.clear()

    def add_result(
        self,
        rank: int,
        sequence: str,
        position: str,
        score: float,
        passed: bool,
        cluster_size: int = 1,
        risk_level: str = "low",
        risk_score: float = 0,
        top_targets: str = "",
        validation_direction: str = "",
        result_index: int = None,
    ):
        """添加一行结果"""
        sorting_enabled = self.table.isSortingEnabled()
        self.table.setSortingEnabled(False)
        row = self.table.rowCount()
        self.table.insertRow(row)

        rank_item = QTableWidgetItem(str(rank))
        rank_item.setData(Qt.ItemDataRole.UserRole, row if result_index is None else result_index)
        self.table.setItem(row, 0, rank_item)
        self.table.setItem(row, 1, QTableWidgetItem(sequence))
        self.table.setItem(row, 2, QTableWidgetItem(position))
        self.table.setItem(row, 3, QTableWidgetItem(f"{score:.2f}"))
        self.table.setItem(row, 4, QTableWidgetItem(str(cluster_size)))
        self.table.setItem(row, 5, QTableWidgetItem(risk_level.title()))
        self.table.setItem(row, 6, QTableWidgetItem(f"{risk_score:.1f}"))
        self.table.setItem(row, 7, QTableWidgetItem(top_targets))
        self.table.setItem(row, 8, QTableWidgetItem(validation_direction))
        self.table.setItem(row, 9, QTableWidgetItem("Yes" if passed else "No"))
        self.table.setSortingEnabled(sorting_enabled)

    def load_results(self, results: list):
        """批量加载结果"""
        sorting_enabled = self.table.isSortingEnabled()
        self.table.setSortingEnabled(False)
        self.clear_results()
        self.results = list(results)
        for result_index, r in enumerate(results):
            off_target = r.get("off_target") or {}
            if not off_target and ("risk_level" in r or "risk_score" in r):
                off_target = {
                    "risk_level": r.get("risk_level", "low"),
                    "risk_score": r.get("risk_score", 0),
                    "top_targets": [
                        {"target_id": item.strip()}
                        for item in str(r.get("top_risk_targets", "")).split(";")
                        if item.strip()
                    ],
                    "validation_direction": r.get("validation_direction", ""),
                }
            top_targets = "; ".join(
                target.get("target_id", "")
                for target in off_target.get("top_targets", [])[:3]
            )
            self.add_result(
                rank=r.get("rank", 0),
                sequence=r.get("sequence", r.get("candidate_seq", "")),
                position=r.get("position", f"{r.get('position_start', '')}-{r.get('position_end', '')}"),
                score=r.get("consensus_score", 0.0),
                passed=r.get("passed", bool(r.get("passed_filters", False))),
                cluster_size=int(r.get("cluster_size", 1) or 1),
                risk_level=off_target.get("risk_level", "low"),
                risk_score=off_target.get("risk_score", 0),
                top_targets=top_targets,
                validation_direction=off_target.get("validation_direction", ""),
                result_index=result_index,
            )
        self.table.setSortingEnabled(sorting_enabled)

    def _show_selected_detail(self):
        """Render details for the selected recommendation."""
        selected = self.table.selectionModel().selectedRows()
        if not selected:
            return
        row = selected[0].row()
        index_item = self.table.item(row, 0)
        result_index = index_item.data(Qt.ItemDataRole.UserRole) if index_item else row
        if result_index is None:
            result_index = row
        result_index = int(result_index)
        if result_index < 0 or result_index >= len(self.results):
            return
        result = self.results[result_index]
        explanation = result.get("explanation") or {}
        validation_hits = result.get("validation_hits") or []
        primers = result.get("primers") or {}

        lines = []
        if explanation:
            lines.append(explanation.get("summary", ""))
            for section in ["efficacy_notes", "risk_notes", "method_notes", "validation_notes"]:
                for item in explanation.get(section, [])[:5]:
                    lines.append(f"- {item}")
        else:
            off_target = result.get("off_target") or {}
            lines.append(f"Risk: {off_target.get('risk_level', 'low')} score {off_target.get('risk_score', 0)}")

        if validation_hits:
            lines.append("")
            lines.append("Off-target validation:")
            for hit in validation_hits[:5]:
                lines.append(
                    f"- {hit.get('target_id', '')}: {hit.get('match_type', '')}, "
                    f"longest {hit.get('longest_contiguous_match', '')}, "
                    f"{hit.get('validation_action', '')}"
                )

        if result.get("region_map"):
            lines.append("")
            lines.append(result["region_map"])

        sgrna = result.get("sgrna") or {}
        if sgrna:
            off_summary = ((result.get("off_target") or {}).get("summary") or {})
            mismatch_counts = off_summary.get("mismatch_counts") or {}
            pot_counts = off_summary.get("pot_mismatch_counts") or {}
            pam = sgrna.get("pam", "")
            genomic_pam = sgrna.get("genomic_pam", pam)
            strand = sgrna.get("strand", "")
            cut_site = sgrna.get("cut_site", "")
            lines.append("")
            lines.append("sgRNA target:")
            lines.append(
                f"- PAM {pam}"
                f"{f' (genomic {genomic_pam})' if genomic_pam and genomic_pam != pam else ''}, "
                f"strand {strand}, cut site {cut_site}"
            )
            if sgrna.get("gc_percent") != "":
                lines.append(f"- Spacer GC {sgrna.get('gc_percent')}%, on-target score {sgrna.get('on_target_score', '')}")
            if mismatch_counts:
                lines.append(
                    "- OT counts 0M-5M: "
                    + ", ".join(f"{bucket}={mismatch_counts.get(bucket, 0)}" for bucket in ["0M", "1M", "2M", "3M", "4M", "5M"])
                )
            if pot_counts:
                lines.append(
                    "- POT seed12 counts 0M-5M: "
                    + ", ".join(f"{bucket}={pot_counts.get(bucket, 0)}" for bucket in ["0M", "1M", "2M", "3M", "4M", "5M"])
                )

        if primers:
            primer_rows = []
            for key in ["forward_primer", "reverse_primer", "t7_forward_primer", "t7_reverse_primer"]:
                primer = primers.get(key) or {}
                if primer.get("sequence"):
                    primer_rows.append((key, primer["sequence"]))
            cloning = primers.get("sgrna_cloning_oligos") or {}
            for key in ["px330_forward_oligo", "px330_reverse_oligo", "custom_forward_oligo", "custom_reverse_oligo"]:
                oligo = cloning.get(key) or {}
                if oligo.get("sequence"):
                    primer_rows.append((key, oligo["sequence"]))
            genotyping = primers.get("genotyping_primers") or {}
            for key in ["forward_primer", "reverse_primer"]:
                primer = genotyping.get(key) or {}
                if primer.get("sequence"):
                    primer_rows.append((f"genotyping_{key}", primer["sequence"]))

        if primers and primer_rows:
            lines.append("")
            lines.append("Primers / oligos:")
            for key, sequence in primer_rows:
                lines.append(f"- {key}: {sequence}")

        self.detail_edit.setPlainText("\n".join(line for line in lines if line is not None))
