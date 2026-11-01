# -*- coding: utf-8 -*-
"""Alternative filler-efficiency formulations: strength-only (primary), porosity-coupled, and (fit only) size-dependent.   python run_alt_cal.py -> shell_cal/alternatives.json"""
import json

import numpy as np

import calibrated as m

q = lambda a: [float(np.percentile(a, p)) for p in (10, 50, 90)]
V = json.load(open("shell_out/validation.json")); VS = json.load(open("shell_out/validation_size.json"))
out = {"fit": {"strength_only": {s: {k: v for k, v in V["per_study_kf"][s].items()} for s in V["per_study_kf"]}, "size_dependent": {s: {k: v for k, v in VS["per_study_kf"][s].items()} for s in VS["per_study_kf"]}}, "claims": {}}
for name, couple in (("strength_only", False), ("porosity_coupled", True)):
    m.KF_COUPLE_POROSITY = couple
    out["claims"][name] = {}
    for case in "ABC":
        d = m.sample(case, 3000, 71); opc = m.evaluate(case, d, dict(path="OPC")); row = {}
        for r, mt in ((0.25, False), (0.25, True), (0.5, True)):
            s = m.summarize(m.evaluate(case, d, dict(path="P1", r=r, d_um=30.0, match=mt), ref=opc), opc)
            row["r%d%s" % (round(100 * r), "m" if mt else "u")] = dict(dE_pct=s["dE_pct"], dL=s["dL"], p_adv_mat=s["p_adv_mat"], p_L=s["p_L"], fc28=s["fc28"], D28_ratio=None)
        sh = m.evaluate(case, d, dict(path="P1", r=0.25, d_um=30.0), ref=opc); row["D28_ratio_r25u"] = q(sh["D28"] / opc["D28"])
        out["claims"][name][case] = row
m.KF_COUPLE_POROSITY = False
json.dump(out, open("shell_cal/alternatives.json", "w"), indent=1); print("ok")
