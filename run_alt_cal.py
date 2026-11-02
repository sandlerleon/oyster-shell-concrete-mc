# -*- coding: utf-8 -*-
"""Alternative filler formulations, transport-model scenarios and leave-one-study-out service-life conclusions (primary draws of the study: seed 21, n = 3000).
    python run_alt_cal.py -> shell_cal/alternatives.json, shell_cal/loso_claims.json"""
import json

import numpy as np

import calibrated as m

q = lambda a: [float(np.percentile(a, p)) for p in (10, 50, 90)]
V = json.load(open("shell_out/validation.json")); VS = json.load(open("shell_out/validation_size.json"))
CFG = (("r25u", 0.25, False), ("r25m", 0.25, True), ("r50m", 0.5, True))


def claims(couple, chi=None, n=3000, seed=21):
    m.KF_COUPLE_POROSITY = couple; out = {}
    for case in "ABC":
        d = m.sample(case, n, seed)
        if chi is not None:
            d["chi_D"] = np.full(n, chi)
        opc = m.evaluate(case, d, dict(path="OPC")); row = {}
        for k, r, mt in CFG:
            s = m.summarize(m.evaluate(case, d, dict(path="P1", r=r, d_um=30.0, match=mt), ref=opc), opc)
            row[k] = dict(dE_pct=s["dE_pct"], dL=s["dL"], p_adv_mat=s["p_adv_mat"], p_L=s["p_L"], fc28=s["fc28"])
        sh = m.evaluate(case, d, dict(path="P1", r=0.25, d_um=30.0), ref=opc); row["D28_ratio_r25u"] = q(sh["D28"] / opc["D28"]); out[case] = row
    m.KF_COUPLE_POROSITY = False
    return out


out = {"fit": {"strength_only": V["per_study_kf"], "size_dependent": VS["per_study_kf"]}, "claims": {"strength_only": claims(False), "porosity_coupled": claims(True)}, "transport": {}}
for chi in (0.0, 0.5, 1.0):
    for couple in (False, True):
        out["transport"]["chi%.1f|%s" % (chi, "coupled" if couple else "strength_only")] = claims(couple, chi)
json.dump(out, open("shell_cal/alternatives.json", "w"), indent=1)

# leave-one-study-out: pooled distribution from the other three studies, full service-life conclusions
ks = V["kf_study_estimates"]; names = list(ks)
loso = {}
orig = m.COMMON["kf"]
for held in names + ["none"]:
    keep = np.array([ks[s] for s in names if s != held]); n = len(keep)
    mu, sd = keep.mean(), keep.std(ddof=1); scale, df = sd * np.sqrt(1 + 1 / n), n - 1
    m.COMMON["kf"] = ("tre", float(mu), float(scale), df, -1.0, 1.5)
    loso[held] = dict(kf=dict(mu=float(mu), scale=float(scale), df=df, n=n), claims=claims(False))
m.COMMON["kf"] = orig
json.dump(loso, open("shell_cal/loso_claims.json", "w"), indent=1)
print("ok")
