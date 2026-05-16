# dsRNA-Forge UI Handoff

## UI Source Files
- `dsforge/gui/main_window.py`: main layout, tab wiring, start/cancel/export/project actions.
- `dsforge/gui/transcript_panel.py`: transcriptome loading, saved transcriptomes, target source, pasted/uploaded targets, extra backgrounds.
- `dsforge/gui/config_panel.py`: design type, confidence preset, advanced parameters, start/cancel buttons.
- `dsforge/gui/results_panel.py`: result table, detail panel, export/project buttons.
- `dsforge/gui/progress_panel.py`: status/progress/log/stats.
- `dsforge/gui/history_panel.py`: historical task table and task loading.
- `dsforge/gui/cache_dialog.py`: saved transcriptome cache management.
- `dsforge/gui/workers.py`: GUI background worker and cancellation.

## Current Main Layout
- Window title: `dsRNA-Forge v0.1.0`
- Minimum size: `1400 x 900`
- Left vertical splitter:
  - `1. Load transcriptome and choose target`
  - `2. Choose design type`
- Right vertical splitter:
  - tabs: `Results`, `History`
  - progress panel
- Initial status bar: `Ready — Load a transcriptome to begin`

## Transcript Panel
- Group title: `1. Load transcriptome and choose target`
- Saved transcriptomes row:
  - Label: `Saved transcriptomes:`
  - Button: `Load Saved`
  - Button: `Manage Cache`
- Transcriptome file row:
  - Placeholder: `Choose transcriptome FASTA...`
  - Button: `Browse...`
  - Load button: `Load & Index`
  - Clear button: `Clear`
- Target source combo:
  - `Search in transcriptome`
  - `Paste target sequence`
  - `Upload target FASTA`
- Search target controls:
  - Placeholder: `Type a gene/transcript ID to filter the target list...`
  - Label: `Target:`
- Paste target controls:
  - Placeholder: `Paste a target sequence here. DNA/RNA and FASTA text are both OK.`
- Upload target controls:
  - Placeholder: `Choose a single-target FASTA...`
  - Button: `Upload Target FASTA...`
- Extra backgrounds:
  - Label format: `Extra backgrounds: N`
  - Button: `Add Background FASTA...`
  - Button: `Clear Backgrounds`
- Stats default: `No transcriptome loaded`
- Log placeholder: `Loading messages will appear here.`

## Config Panel
- Group title: `2. Choose design type`
- Design type combo:
  - `siRNA (21 nt)`
  - `DsiRNA (27 nt)`
  - `Long dsRNA for RNAi (200-500 bp)`
  - `sgRNA for SpCas9 (20 nt + NGG PAM)`
- Design confidence combo:
  - `Strict - lowest off-target risk`
  - `Balanced - recommended`
  - `Relaxed - rescue difficult targets`
- Advanced toggle: `Show advanced settings`
- Advanced group: `Advanced Settings`
  - `Length Range`: `Min:`, `Max:`
  - `GC Content (%)`: `Min:`, `Max:`
  - `Scoring Rules`: `Consensus (required)`, `Reynolds`, `Ui-Tei`, `Amarzguioui`, `Hsieh`, `Jagla`
  - `CPU cores:`
- Buttons:
  - `Start Design`
  - `Cancel`

## Results Panel
- Group title: `Results`
- Table columns:
  - `Rank`
  - `Sequence`
  - `Position`
  - `Consensus Score`
  - `Cluster Size`
  - `Risk`
  - `Risk Score`
  - `Top Risk Targets`
  - `Validation Direction`
  - `Pass`
- Detail placeholder:
  - `Select a recommendation to inspect ranking reasons, off-target validation, region map and primers/oligos.`
- Detail sections currently rendered:
  - explanation summary
  - efficacy/risk/method/validation notes
  - `Off-target validation:`
  - ASCII region map
  - `sgRNA target:`
  - sgRNA OT counts: `OT counts 0M-5M: 0M=..., 1M=..., ...`
  - sgRNA POT counts: `POT seed12 counts 0M-5M: 0M=..., 1M=..., ...`
  - `Primers / oligos:`
- Buttons:
  - `Export CSV`
  - `Export Excel`
  - `Export FASTA`
  - `Export Report`
  - `Export Primers`
  - `Save Project`
  - `Open Project`

## History And Progress
- History panel source file: `dsforge/gui/history_panel.py`
- Progress panel source file: `dsforge/gui/progress_panel.py`
- Runtime phases emitted by design tasks include candidate generation, off-target index build/reuse, scoring, ranking, RNAup refinement, saving, and table loading.

## Export Entry Points
- CSV: `ResultExporter.export_csv`
- Excel: `ResultExporter.export_excel`
- FASTA: `ResultExporter.export_fasta`
- Validation report: `ResultExporter.export_validation_report`
- Primer/oligo order CSV: `ResultExporter.export_primer_order_csv`

## Phase22 dsRNA Off-target Report Fields
- `result["off_target"]["top_targets"]` is the source of potential RNAi off-target genes/transcripts for siRNA, DsiRNA and long dsRNA results.
- CSV field names:
  - `top_risk_targets`: compact `target_id:risk_score` list.
  - `top_offtarget_genes`: user-facing list with target ID, risk score and top risk reasons.
- Excel columns:
  - `Top Risk Targets`
  - `Top Off-target Genes`
- Validation report:
  - `Recommendations` sheet now includes `Top Off-target Genes` so the first sheet explicitly tells the user which potential off-target genes/transcripts to validate.
  - `OffTarget Validation` sheet includes `Reasons` next to each target hit.

## Phase20 sgRNA Fields
- Design PAM remains SpCas9 canonical `NGG`; off-target search now checks SpCas9 `NRG` sites, including `NAG`.
- Off-target summary keys under `result["off_target"]["summary"]`:
  - `mismatch_counts`: `0M`, `1M`, `2M`, `3M`, `4M`, `5M`
  - `pot_mismatch_counts`: same buckets, only hits with 12 nt PAM-proximal seed match
  - `total_ot`
  - `total_pot`
  - `nrg_pam_hits`: all NRG off-target hits counted in the sgRNA risk scan
  - `nag_pam_hits`: alternative NAG PAM hits within those NRG hits
  - `sgrnacas9_risk_evaluation`
- CSV field names:
  - `sgrna_ot_0M` ... `sgrna_ot_5M`
  - `sgrna_pot_0M` ... `sgrna_pot_5M`
  - `sgrna_total_ot`
  - `sgrna_total_pot`
  - `sgrna_nrg_pam_hits`
  - `sgrna_nag_pam_hits`
  - `sgrna_risk_evaluation`
- Cloning oligo keys under `result["primers"]["sgrna_cloning_oligos"]`:
  - existing: `px330_forward_oligo`, `px330_reverse_oligo`
  - new: `custom_forward_oligo`, `custom_reverse_oligo`
  - warning list: `restriction_site_warnings`

## Minimal-Change UI Notes
- UI code is mostly isolated under `dsforge/gui/`.
- Keep signal names stable when possible:
  - `transcriptome_loaded`
  - `transcriptome_cleared`
  - `manage_cache_requested`
  - worker signals: `progress`, `result`, `error`, `finished`
- Keep `ResultsPanel.load_results(results)` and `ConfigPanel.get_config()` return shapes stable unless also updating `MainWindow`.
- Algorithm, database, export, and runtime packaging logic live outside `dsforge/gui/`; UI polish should usually not touch them.
