# -*- coding: utf-8 -*-
"""Prior-predictive relative strength of the ORIGINAL model (kf = 0, inert filler) for the three published mixes; python kf0_prediction.py -> shell_out/kf0_pred.json"""
import json

import numpy as np

import ShellConcreteMC_v3 as m
MIXES = {"S0F15": dict(r=50 / 350, d_um=90.0), "S15F15": dict(r=95 / 350, d_um=64.0), "S30F15": dict(r=140 / 350, d_um=55.0)}
OBS = {"S0F15": 1.00, "S15F15": (1 - 0.07) ** 2, "S30F15": 0.75}
TOL = 0.07


def rel_strength(d, r, d_um):
    n = len(d["B"]); B = np.full(n, 350.0); w = np.full(n, 160.0); T = np.full(n, 20.0)
    return m.fc28(d, B, w, r, d_um, T) / m.fc28(d, B, w, 0.0, d_um, T)


out = {}
for case in "ABC":
    d = m.sample(case, 100000, seed=101)
    d["kf"] = np.zeros(100000)
    rel = {k: rel_strength(d, v["r"], v["d_um"]) for k, v in MIXES.items()}
    acc = np.ones(100000, bool)
    for k in ("S0F15", "S30F15"):
        acc &= np.abs(rel[k] - OBS[k]) <= TOL
    out[case] = dict(pred={k: [float(np.percentile(x, p)) for p in (10, 50, 90)] for k, x in rel.items()},
                     in_tol={k: float(np.mean(np.abs(rel[k] - OBS[k]) <= TOL)) for k in rel}, joint_cal_accept=int(acc.sum()))
json.dump(out, open("shell_out/kf0_pred.json", "w"), indent=1)
print(json.dumps(out["A"], indent=1))
