"""
参数配置面板
- 设计模式选择
- 规则开关
- 长度/GC 范围
- 脱靶/热力学开关
- 核数选择
"""

import os
from PyQt6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QComboBox,
    QSpinBox,
    QCheckBox,
    QGroupBox,
    QPushButton,
    QGridLayout,
)
from PyQt6.QtCore import Qt


class ConfigPanel(QGroupBox):
    """参数配置面板"""

    def __init__(self):
        super().__init__("2. Choose design type")
        self._setup_ui()
        self._load_defaults()

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(10)

        # === 设计模式 ===
        mode_layout = QHBoxLayout()
        mode_layout.addWidget(QLabel("Design type:"))
        self.mode_combo = QComboBox()
        self.mode_combo.addItems([
            "siRNA (21 nt)",
            "DsiRNA (27 nt)",
            "Long dsRNA for RNAi (200-500 bp)",
            "sgRNA for SpCas9 (20 nt + NGG PAM)",
        ])
        self.mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        mode_layout.addWidget(self.mode_combo)
        layout.addLayout(mode_layout)

        preset_layout = QHBoxLayout()
        preset_layout.addWidget(QLabel("Design confidence:"))
        self.preset_combo = QComboBox()
        self.preset_combo.addItems([
            "Strict - lowest off-target risk",
            "Balanced - recommended",
            "Relaxed - rescue difficult targets",
        ])
        self.preset_combo.setCurrentIndex(1)
        self.preset_combo.currentIndexChanged.connect(self._on_preset_changed)
        preset_layout.addWidget(self.preset_combo)
        layout.addLayout(preset_layout)

        self.advanced_toggle = QCheckBox("Show advanced settings")
        self.advanced_toggle.toggled.connect(self._on_advanced_toggled)
        layout.addWidget(self.advanced_toggle)

        self.advanced_group = QGroupBox("Advanced Settings")
        advanced_layout = QVBoxLayout(self.advanced_group)

        # === 长度范围 ===
        len_group = QGroupBox("Length Range")
        len_layout = QGridLayout(len_group)

        len_layout.addWidget(QLabel("Min:"), 0, 0)
        self.min_len_spin = QSpinBox()
        self.min_len_spin.setRange(15, 2000)
        len_layout.addWidget(self.min_len_spin, 0, 1)

        len_layout.addWidget(QLabel("Max:"), 1, 0)
        self.max_len_spin = QSpinBox()
        self.max_len_spin.setRange(15, 2000)
        len_layout.addWidget(self.max_len_spin, 1, 1)

        advanced_layout.addWidget(len_group)

        # === GC 含量 ===
        gc_group = QGroupBox("GC Content (%)")
        gc_layout = QGridLayout(gc_group)

        gc_layout.addWidget(QLabel("Min:"), 0, 0)
        self.min_gc_spin = QSpinBox()
        self.min_gc_spin.setRange(0, 100)
        gc_layout.addWidget(self.min_gc_spin, 0, 1)

        gc_layout.addWidget(QLabel("Max:"), 1, 0)
        self.max_gc_spin = QSpinBox()
        self.max_gc_spin.setRange(0, 100)
        gc_layout.addWidget(self.max_gc_spin, 1, 1)

        advanced_layout.addWidget(gc_group)

        # === 规则开关 ===
        rules_group = QGroupBox("Scoring Rules")
        rules_layout = QVBoxLayout(rules_group)

        self.rule_consensus = QCheckBox("Consensus (required)")
        self.rule_consensus.setChecked(True)
        self.rule_consensus.setEnabled(False)
        rules_layout.addWidget(self.rule_consensus)

        self.rule_reynolds = QCheckBox("Reynolds")
        rules_layout.addWidget(self.rule_reynolds)

        self.rule_ui_tei = QCheckBox("Ui-Tei")
        rules_layout.addWidget(self.rule_ui_tei)

        self.rule_amarzguioui = QCheckBox("Amarzguioui")
        rules_layout.addWidget(self.rule_amarzguioui)

        self.rule_hsieh = QCheckBox("Hsieh")
        rules_layout.addWidget(self.rule_hsieh)

        self.rule_jagla = QCheckBox("Jagla")
        rules_layout.addWidget(self.rule_jagla)

        advanced_layout.addWidget(rules_group)

        # === 并行设置 ===
        para_layout = QHBoxLayout()
        para_layout.addWidget(QLabel("CPU cores:"))
        self.cores_spin = QSpinBox()
        self.cores_spin.setRange(1, os.cpu_count() or 4)
        self.cores_spin.setValue(max(1, (os.cpu_count() or 4) - 1))
        para_layout.addWidget(self.cores_spin)
        advanced_layout.addLayout(para_layout)
        layout.addWidget(self.advanced_group)
        self.advanced_group.setVisible(False)

        # === 操作按钮 ===
        btn_layout = QHBoxLayout()
        self.start_btn = QPushButton("Start Design")
        self.start_btn.setStyleSheet("font-weight: bold;")
        self.start_btn.setEnabled(False)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)

        btn_layout.addWidget(self.start_btn)
        btn_layout.addWidget(self.cancel_btn)
        layout.addLayout(btn_layout)

        layout.addStretch()

    def _load_defaults(self):
        """加载默认配置"""
        self._on_mode_changed(self.mode_combo.currentIndex())
        self.min_gc_spin.setValue(30)
        self.max_gc_spin.setValue(60)
        self.rule_reynolds.setChecked(True)
        self.rule_ui_tei.setChecked(True)

    def _on_mode_changed(self, index: int):
        """模式切换时更新长度范围"""
        ranges = [
            (21, 21),      # siRNA
            (27, 27),      # DsiRNA
            (200, 500),    # Long dsRNA
            (20, 20),      # sgRNA spacer
        ]
        min_len, max_len = ranges[index]
        self.min_len_spin.setValue(min_len)
        self.max_len_spin.setValue(max_len)

    def _on_preset_changed(self, index: int):
        """Adjust visible simple defaults when the user changes preset."""
        if index == 0:
            self.min_gc_spin.setValue(35)
            self.max_gc_spin.setValue(50)
        elif index == 2:
            self.min_gc_spin.setValue(20)
            self.max_gc_spin.setValue(65)
        else:
            self.min_gc_spin.setValue(30)
            self.max_gc_spin.setValue(60)

    def _on_advanced_toggled(self, checked: bool):
        """Show or hide expert parameters."""
        self.advanced_group.setVisible(checked)

    def get_config(self) -> dict:
        """获取当前配置"""
        modes = ["siRNA", "DsiRNA", "long_dsRNA", "sgRNA"]
        presets = ["strict", "balanced", "relaxed"]
        return {
            "mode": modes[self.mode_combo.currentIndex()],
            "preset": presets[self.preset_combo.currentIndex()],
            "length": {
                "min": self.min_len_spin.value(),
                "max": self.max_len_spin.value(),
            },
            "gc": {
                "min": self.min_gc_spin.value(),
                "max": self.max_gc_spin.value(),
            },
            "rules": {
                "consensus": True,
                "reynolds": self.rule_reynolds.isChecked(),
                "ui_tei": self.rule_ui_tei.isChecked(),
                "amarzguioui": self.rule_amarzguioui.isChecked(),
                "hsieh": self.rule_hsieh.isChecked(),
                "jagla": self.rule_jagla.isChecked(),
            },
            "cores": self.cores_spin.value(),
        }

    def apply_config(self, config: dict):
        """Restore UI controls from a saved DesignConfig-like dict."""
        if not isinstance(config, dict):
            return

        modes = ["siRNA", "DsiRNA", "long_dsRNA", "sgRNA"]
        mode = config.get("mode")
        if mode in modes:
            self.mode_combo.setCurrentIndex(modes.index(mode))

        presets = ["strict", "balanced", "relaxed"]
        preset = config.get("preset", "balanced")
        if preset in presets:
            self.preset_combo.setCurrentIndex(presets.index(preset))

        length = config.get("length") or {}
        self.min_len_spin.setValue(int(config.get("length_min", length.get("min", self.min_len_spin.value()))))
        self.max_len_spin.setValue(int(config.get("length_max", length.get("max", self.max_len_spin.value()))))

        gc = config.get("gc") or {}
        self.min_gc_spin.setValue(int(config.get("gc_min", gc.get("min", self.min_gc_spin.value()))))
        self.max_gc_spin.setValue(int(config.get("gc_max", gc.get("max", self.max_gc_spin.value()))))

        enabled_rules = config.get("enabled_rules")
        rules = config.get("rules") or {}
        if isinstance(enabled_rules, list):
            rules = {name: name in enabled_rules for name in ["reynolds", "ui_tei", "amarzguioui", "hsieh", "jagla"]}
        self.rule_reynolds.setChecked(bool(rules.get("reynolds", self.rule_reynolds.isChecked())))
        self.rule_ui_tei.setChecked(bool(rules.get("ui_tei", self.rule_ui_tei.isChecked())))
        self.rule_amarzguioui.setChecked(bool(rules.get("amarzguioui", self.rule_amarzguioui.isChecked())))
        self.rule_hsieh.setChecked(bool(rules.get("hsieh", self.rule_hsieh.isChecked())))
        self.rule_jagla.setChecked(bool(rules.get("jagla", self.rule_jagla.isChecked())))

        cores = int(config.get("n_cores", config.get("cores", self.cores_spin.value())) or 1)
        self.cores_spin.setValue(max(self.cores_spin.minimum(), min(self.cores_spin.maximum(), cores)))
