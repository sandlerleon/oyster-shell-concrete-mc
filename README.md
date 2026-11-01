# Oyster-shell concrete Monte Carlo

Leon Sandler, Independent Researcher (Northbrook, Illinois, USA) - sandler.leon@gmail.com - ORCID 0009-0007-4584-808X

Code, tests, results and manuscript for *"A Coupled Chemo-Mechanical and Life-Cycle Framework for Oyster-Shell-Derived Cementitious Materials in Aggressive Infrastructure Environments"*, prepared for **Construction and Building Materials** (Elsevier).

Archived on Zenodo: software doi:[10.5281/zenodo.23867048](https://doi.org/10.5281/zenodo.23867048) (all versions), manuscript preprint doi:[10.5281/zenodo.23867044](https://doi.org/10.5281/zenodo.23867044) (all versions).
Companion paper (allocation of shell between reef restoration, concrete and disposal; Journal of Cleaner Production): [https://github.com/sandlerleon/oyster-shell-allocation-mc](https://github.com/sandlerleon/oyster-shell-allocation-mc), doi:[10.5281/zenodo.23867050](https://doi.org/10.5281/zenodo.23867050).

**This is an in-silico study.** No new experiments are reported. Every probability is conditional on the stated priors and on the model structure. Only the strength model has been tested against published data (one scallop-shell dataset, three mixes); chloride, acid and freeze-thaw transport remain uncalibrated.

## What it does

A coupled model links hydration kinetics and porosity, chloride, acid and leaching transport, damage accumulation (fatigue, freeze-thaw, shrinkage, thermal, abrasion) and cradle-to-gate carbon for two pathways (P1 shell powder replacing cement; P2 shell CaCO3 replacing limestone in kiln feed) in three cases (pavement, dam, marine culvert). It is propagated by Latin hypercube Monte Carlo with paired shell/OPC draws, Sobol sensitivity and a falsification region.

**Out-of-sample test (v3).** `calibrate_abbas.py` tests the strength model against Abbas and Kudukkan (2025, Constr. Mater. 5(4):82). The original inert-filler model is rejected (`kf0_prediction.py`); a filler-efficiency parameter `kf` (v3 extension) is calibrated by rejection ABC on two mixes with the third held out, and the main claims are re-evaluated under the posterior.

## Layout

| Path | Content |
|---|---|
| `ShellConcreteMC_v2.py` | original model (unchanged) with `--verify` checks |
| `ShellConcreteMC_v3.py` | v2 plus the optional filler efficiency `kf` (`kf = 0` reproduces v2 exactly) |
| `run_stage1.py`, `run_sobol.py`, `run_extras.py`, `run_figs.py` | stages that regenerate the results and figures of the original analysis |
| `calibrate_abbas.py`, `kf0_prediction.py`, `make_fig_calibration.py` | out-of-sample test, calibration and its figures |
| `shell_out/` | all result files (JSON) and figures |
| `tests/test_p1.py` | automated checks |
| `manuscript/` | manuscript (docx, PDF), supplementary information, cover letter, highlights, builder script |

## Reproduce

```
pip install -r requirements.txt
python ShellConcreteMC_v2.py --verify
python run_stage1.py && python run_sobol.py A 2048 && python run_sobol.py B 2048 && python run_sobol.py C 2048 && python run_extras.py && python run_figs.py
python calibrate_abbas.py && python kf0_prediction.py && python make_fig_calibration.py shell_out
python -m pytest tests
```

The original analysis was re-run from a clean directory and reproduced the published result files to numerical precision (maximum absolute difference about 1e-13).

## Limitations

Literature-style priors; reference lives anchored to plausible design lives (results are relative to the OPC reference); one calibration dataset for scallop (not oyster) shell with strengths partly read from text; transport and damage uncalibrated; plain-OPC reference; generic emission factors.

## Licence

MIT (see `LICENSE`).
