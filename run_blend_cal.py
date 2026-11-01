# -*- coding: utf-8 -*-
"""Shell concrete (P1) against a plain-OPC reference and against a blended-cement reference (50 % slag) under the calibrated filler efficiency.
    python run_blend_cal.py  ->  shell_cal/blend.json"""
import json

import numpy as np

import calibrated as m

q = lambda a: [float(np.percentile(a, p)) for p in (10, 50, 90)]
out = {"fs": 0.5, "EF_slag": 0.07, "cases": {}}
for case in "ABC":
    d = m.sample(case, 3000, 61)
    opc = m.evaluate(case, d, dict(path="OPC"))
    row = {"opc_L": q(opc["L"]), "opc_E": q(opc["E_gate"])}
    for tag, cfg in (("blend50", dict(path="BLEND", fs=0.5)), ("blend30", dict(path="BLEND", fs=0.3)), ("blend50_slag12", dict(path="BLEND", fs=0.5, EF_slag=0.12))):
        bl = m.evaluate(case, d, cfg, ref=opc)
        row[tag] = dict(L=q(bl["L"]), dL_vs_opc=q(bl["L"] - opc["L"]), dE_pct_vs_opc=q(100 * (bl["E_gate"] / opc["E_gate"] - 1)), fc28_ratio=float(np.median(bl["fc28"] / opc["fc28"])))
    for r in (0.25, 0.50):
        for match in (False, True):
            sh = m.evaluate(case, d, dict(path="P1", r=r, d_um=30.0, match=match), ref=opc)
            for tag, cfg in (("blend50", dict(path="BLEND", fs=0.5)), ("blend30", dict(path="BLEND", fs=0.3))):
                bl = m.evaluate(case, d, cfg, ref=opc)
                s = m.summarize(sh, bl)
                row["P1_r%d%s_vs_%s" % (round(100 * r), "m" if match else "u", tag)] = dict(dE_pct=s["dE_pct"], dL=s["dL"], p_adv_mat=s["p_adv_mat"], p_L=s["p_L"], p_E=s["p_E"])
            s0 = m.summarize(sh, opc)
            row["P1_r%d%s_vs_opc" % (round(100 * r), "m" if match else "u")] = dict(dE_pct=s0["dE_pct"], dL=s0["dL"], p_adv_mat=s0["p_adv_mat"], p_L=s0["p_L"], p_E=s0["p_E"])
    out["cases"][case] = row
    print(case, "ok", flush=True)
json.dump(out, open("shell_cal/blend.json", "w"), indent=1)
