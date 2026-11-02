# Oyster-shell concrete Monte Carlo

Leon Sandler, Independent Researcher (Northbrook, Illinois, USA) - sandler.leon@gmail.com - ORCID 0009-0007-4584-808X

Code, tests, results and manuscript for *"A Coupled Chemo-Mechanical and Life-Cycle Framework for Oyster-Shell-Derived Cementitious Materials in Aggressive Infrastructure Environments"*, prepared for **Construction and Building Materials** (Elsevier). Current version: 1.2.3. Version 1.2.3 changes only the manuscript files: the language-edited text was merged selectively (see `manuscript/README_MANUSCRIPT.md`) and the keyword list was reduced to six; no code, result or figure changed.

Archived on Zenodo: software doi:[10.5281/zenodo.23867048](https://doi.org/10.5281/zenodo.23867048) (all versions), manuscript preprint doi:[10.5281/zenodo.23867044](https://doi.org/10.5281/zenodo.23867044) (all versions).
Companion paper (allocation of shell between reef restoration, concrete and disposal; Journal of Cleaner Production): [https://github.com/sandlerleon/oyster-shell-allocation-mc](https://github.com/sandlerleon/oyster-shell-allocation-mc), doi:[10.5281/zenodo.23867050](https://doi.org/10.5281/zenodo.23867050).

**This is an in-silico study.** No new experiments are reported. Every probability is conditional on the stated priors and on the model structure. Only the strength model has been tested against published data; chloride, acid and freeze-thaw transport remain uncalibrated.

## What it does

A coupled model links hydration kinetics and porosity, chloride, acid and leaching transport, damage accumulation and cradle-to-gate carbon for two pathways (P1 shell powder replacing cement; P2 shell CaCO3 replacing limestone in kiln feed) in three cases (pavement, dam, marine culvert), with paired shell/OPC Monte Carlo draws, Sobol sensitivity and a falsification region.

**Version 1.2.2 changes.** (1) The inert-filler strength model is rejected by scallop-shell data (`calibrate_abbas.py`, `kf0_prediction.py`). (2) An effective strength-efficiency coefficient `k_f` calibrated on that study is tested, with frozen parameters, on 11 independent oyster-shell points from three studies (`validate_independent.py`); it does not outperform dilution-only or inert models, and `k_f` differs between studies. (3) The primary results use a pooled random-effects distribution of `k_f` (`calibrated.py`, `run_*_cal.py`); the original and single-study variants are ablations. (4) A 50 % slag blended-cement reference (`run_blend_cal.py`) and alternative formulations (porosity-coupled and size-dependent `k_f`; `run_alt_cal.py`) are added. (5) `export_joint_draws.py` writes the joint draws used by the companion paper. (6) Version 1.2.0 adds a table of the conditions of every probability, a leave-one-study-out analysis of the service-life conclusions (`run_alt_cal.py`), fixed transport-model scenarios, study-level prediction errors and a table of calibrated, uncalibrated and assumed components.

## Layout

| Path | Content |
|---|---|
| `ShellConcreteMC_v2.py` | original model (unchanged) with `--verify` checks |
| `ShellConcreteMC_v3.py` | v2 plus the optional coefficient `k_f` (zero reproduces v2), a blended-cement reference path and alternative `k_f` switches |
| `calibrated.py` | v3 with `k_f` ~ pooled predictive distribution (`KF_MODE=zero` for the original model) |
| `calibrate_abbas.py`, `kf0_prediction.py`, `validate_independent.py` | calibration, original-model prediction, independent validation and per-study `k_f` |
| `run_stage1*.py`, `run_sobol*.py`, `run_extras*.py`, `run_figs*.py` | stages of the original (`shell_out/`) and primary (`*_cal`, `shell_cal/`) analyses |
| `run_blend_cal.py`, `run_alt_cal.py`, `make_fig_*.py`, `export_joint_draws.py` | blended-cement benchmark, alternatives, figures, export of joint draws |
| `shell_out/`, `shell_cal/` | result files (JSON) and figures |
| `build_manuscript.py` | regenerates the manuscript from the result files and compares it with the shipped one |
| `tests/test_p1.py` | automated checks |
| `manuscript/` | manuscript (docx, PDF), supplementary information, cover letter, highlights, builder script |

## Reproduce

```
pip install -r requirements.txt
python ShellConcreteMC_v2.py --verify
python calibrate_abbas.py && python kf0_prediction.py && python validate_independent.py
python run_stage1_cal.py && python run_sobol_cal.py A 2048 && python run_sobol_cal.py B 2048 && python run_sobol_cal.py C 2048 && python run_extras_cal.py && python run_figs_cal.py
python run_blend_cal.py && python run_alt_cal.py && python make_fig_calibration.py shell_out && python make_fig_validation.py shell_cal
python build_manuscript.py --compare   # rebuilds the manuscript (needs Node.js) and compares every paragraph and table cell with the shipped one
python -m pytest tests
```

The original analysis was re-run from a clean directory and reproduced the original result files to numerical precision (maximum absolute difference about 1e-13). The environment is pinned in `requirements.txt` (Python 3.11); the random seeds are fixed in the scripts (calibration 101/555, validation 2024/99, studies 21, scans 31, Sobol 11, blended benchmark 61, alternatives and leave-one-study-out 21, joint draws 500-530) and are listed in the manuscript appendix; `REPRODUCTION_CHECK.md` records a clean-directory re-run of the scripts against the shipped result files.

## Limitations

Literature-style priors; reference lives anchored to plausible design lives (results are relative to the reference); calibration data are one scallop-shell study and three independent oyster-shell studies, one of which is a secondary transcription and one with coarse particle fractions; the pooled coefficient describes between-study heterogeneity and not a material property; transport and damage are uncalibrated; the blended-cement reference is a screening treatment; generic emission factors.

## Licence

MIT (see `LICENSE`).
