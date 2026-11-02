# -*- coding: utf-8 -*-
"""Regenerate the manuscript (tables, numbers, figures placement) from the result files of this repository and compare it with the shipped manuscript.

    python build_manuscript.py            build manuscript/_build/Manuscript.docx from shell_cal/, shell_out/ and the figures
    python build_manuscript.py --compare  also compare the text of every paragraph and table cell with manuscript/Oyster_P1_CBM_Manuscript.docx

Requires Node.js (npm packages are installed on first use: docx@9) and python-docx.
"""
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "manuscript", "src"); BUILD = os.path.join(HERE, "manuscript", "_build")
if os.path.exists(BUILD):
    shutil.rmtree(BUILD)
os.makedirs(os.path.join(BUILD, "shell_out")); os.makedirs(os.path.join(BUILD, "shell_orig"))
cal, orig = os.path.join(HERE, "shell_cal"), os.path.join(HERE, "shell_out")
for f in ("results", "extras", "stage1", "sobol_A", "sobol_B", "sobol_C", "blend", "alternatives", "loso_claims"):
    shutil.copy(os.path.join(cal, f + ".json"), os.path.join(BUILD, "shell_out", f + ".json"))
for f in ("params", "calibration", "kf0_pred", "validation", "validation_size"):
    shutil.copy(os.path.join(orig, f + ".json"), os.path.join(BUILD, "shell_out", f + ".json"))
for f in ("results", "extras", "stage1"):                       # original (inert-filler) analysis, used for the ablation
    shutil.copy(os.path.join(orig, f + ".json"), os.path.join(BUILD, "shell_orig", f + ".json"))
shutil.copy(os.path.join(SRC, "dois.json"), os.path.join(BUILD, "shell_out", "dois.json"))
for f in os.listdir(cal):
    if f.startswith("fig") and f.endswith(".png"):
        shutil.copy(os.path.join(cal, f), os.path.join(BUILD, "shell_out", f))
shutil.copy(os.path.join(orig, "fig11_calibration.png"), os.path.join(BUILD, "shell_out", "fig11_calibration.png"))
shutil.copy(os.path.join(SRC, "build_shell_manuscript.js"), os.path.join(BUILD, "build_shell_manuscript.js"))
shell = sys.platform == "win32"
if not os.path.isdir(os.path.join(BUILD, "node_modules", "docx")):
    subprocess.check_call(["npm", "install", "--no-audit", "--no-fund", "docx@9"], cwd=BUILD, shell=shell)
subprocess.check_call(["node", "build_shell_manuscript.js"], cwd=BUILD, shell=shell)
out = os.path.join(BUILD, "Shell_Concrete_Manuscript_MaterialsStructures.docx")
print("built", out)
if "--compare" in sys.argv:
    import docx
    from docx.oxml.ns import qn

    def paras(p):
        return ["".join(x.text or "" for x in e.iter(qn("w:t"))) for e in docx.Document(p).element.body.iter(qn("w:p")) if "".join(x.text or "" for x in e.iter(qn("w:t"))).strip()]
    a, b = paras(out), paras(os.path.join(HERE, "manuscript", "Oyster_P1_CBM_Manuscript.docx"))
    diff = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    print("paragraphs and table cells: built %d, shipped %d, differing %d" % (len(a), len(b), len(diff)))
    sys.exit(0 if len(a) == len(b) and not diff else 1)
