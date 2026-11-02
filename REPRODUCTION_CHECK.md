# Reproduction check (version 1.2.0)

The repository was copied to a clean directory (without the git history) and every script was re-run: `ShellConcreteMC_v2.py --verify`; `calibrate_abbas.py`, `kf0_prediction.py`, `validate_independent.py` (also with `KF_SIZE=1`); the original pipeline (`run_stage1.py`, `run_sobol.py A|B|C 2048`, `run_extras.py`, `run_figs.py`); the primary pipeline (`run_stage1_cal.py`, `run_sobol_cal.py A|B|C 2048`, `run_extras_cal.py`, `run_figs_cal.py`); `run_blend_cal.py`, `run_alt_cal.py` (which also writes `loso_claims.json`) and `export_joint_draws.py`.

Every regenerated JSON file in `shell_out/` and `shell_cal/` was compared with the shipped file: the maximum absolute difference over all numbers is 0 for all 29 files. Seeds are fixed in the scripts and listed in the manuscript appendix. Environment: Python 3.11.9, numpy 2.4.6, scipy 1.17.1, matplotlib 3.11.1 (pinned in `requirements.txt`).

## Manuscript, tables and figures (version 1.2.3; `build_manuscript.py --compare` regenerates the text before the language edit, see manuscript/README_MANUSCRIPT.md)

- `python build_manuscript.py --compare` was run in a clean copy: it stages the result files, runs the builder (`manuscript/src/build_shell_manuscript.js`, Node.js with `docx@9`) and compares the text of every paragraph and table cell with the shipped manuscript: 1,417 of 1,417 identical. The manuscript has 21 tables and 14 figures, all generated from the result files and figure scripts of this repository.
- All 13 figures written to `shell_cal/` by the scripts were regenerated and are byte-identical (MD5) to the shipped files; `make_fig_calibration.py` regenerates `shell_out/fig11_calibration.png` identically.
- The three k_f specifications (original, single-study, pooled) and the seeds are listed in Appendix A and C of the manuscript.
