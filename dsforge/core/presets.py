"""
Experiment-facing design presets.
"""

import copy


def apply_preset_to_config(config, preset: str):
    """Return a config copy adjusted for strict/balanced/relaxed workflows."""
    preset = (preset or "balanced").strip().lower()
    if preset not in {"strict", "balanced", "relaxed"}:
        preset = "balanced"

    cfg = copy.deepcopy(config)
    cfg.preset = preset

    if preset == "strict":
        cfg.gc_min = max(cfg.gc_min, 35.0)
        cfg.gc_max = min(cfg.gc_max, 50.0)
        cfg.off_target_levels = {
            "level_1_16bp": True,
            "level_2_20bp": True,
            "level_3_27bp": True,
            "seed_7nt": True,
        }
        cfg.thermodynamics = {**(cfg.thermodynamics or {}), "rnaup_top_n": 50}
        cfg.cluster_overlap_threshold = min(cfg.cluster_overlap_threshold, 0.75)
    elif preset == "relaxed":
        cfg.gc_min = min(cfg.gc_min, 20.0)
        cfg.gc_max = max(cfg.gc_max, 65.0)
        cfg.off_target_levels = {
            "level_1_16bp": True,
            "level_2_20bp": False,
            "level_3_27bp": False,
            "seed_7nt": True,
        }
        cfg.thermodynamics = {**(cfg.thermodynamics or {}), "rnaup_top_n": 10}
        cfg.cluster_overlap_threshold = max(cfg.cluster_overlap_threshold, 0.85)
    else:
        cfg.gc_min = min(cfg.gc_min, 30.0)
        cfg.gc_max = max(cfg.gc_max, 60.0)
        cfg.off_target_levels = {
            "level_1_16bp": True,
            "level_2_20bp": True,
            "level_3_27bp": False,
            "seed_7nt": True,
        }
        cfg.thermodynamics = {**(cfg.thermodynamics or {}), "rnaup_top_n": 20}

    return cfg
