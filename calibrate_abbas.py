# -*- coding: utf-8 -*-
"""Out-of-sample test and rejection-ABC calibration of the shell-specific strength parameters against the only open quantitative dataset found
(Abbas and Kudukkan 2025, Constr. Mater. 5(4):82, scallop-shell powder, binder 350 kg/m3, water 160 kg/m3, 28-day strength).

    python calibrate_abbas.py   ->  shell_out/calibration.json
"""
import json
import time

import numpy as np

import ShellConcreteMC_v3 as m

T0 = time.time()
# published mixes (kg/m3): control 350 cement; S0F15 300 cement + 50 filler; S15F15 255 + 45 + 50; S30F15 210 + 90 + 50. Binder (cement + shell) is 350 in all mixes; water 160.
MIXES = {"S0F15": dict(r=50 / 350, d_um=90.0), "S15F15": dict(r=95 / 350, d_um=64.0), "S30F15": dict(r=140 / 350, d_um=55.0)}
# observed relative 28-d strength: text ("almost the same as the control", "nearly 25 % reduction") and, for S15F15, derived from the reported moduli of elasticity (-5, -7, -13 %) with E ~ sqrt(fc)
OBS = {"S0F15": 1.00, "S15F15": (1 - 0.07) ** 2, "S30F15": 0.75}
E_RATIO = {"S0F15": 0.95, "S15F15": 0.93, "S30F15": 0.87}
CAL_SET, HOLDOUT = ("S0F15", "S30F15"), "S15F15"
TOL = 0.07


def rel_strength(d, r, d_um):
    n = len(d["B"]); B = np.full(n, 350.0); w = np.full(n, 160.0); T = np.full(n, 20.0)
    return m.fc28(d, B, w, r, d_um, T) / m.fc28(d, B, w, 0.0, d_um, T)


def shell_ef(d, d_um):
    return (d["e_dry"] * d["EF_th"] + d["E_ref"] * (45.0 / d_um) ** 0.5 * d["EF_el"]) / 1000.0 + d["dist_sh"] * d["EF_tr"] / 1000.0 - d["c_land"]


out = {"obs": OBS, "tol": TOL, "mixes": MIXES, "cases": {}}
N = 300000
for case in "ABC":
    d = m.sample(case, N, seed=101)
    d["kf"] = np.random.default_rng(555).uniform(0.0, 1.2, N)      # extension: filler efficiency (k-value), prior U(0, 1.2); kf = 0 is the original model
    rel = {k: rel_strength(d, v["r"], v["d_um"]) for k, v in MIXES.items()}
    prior = {k: [float(np.percentile(x, p)) for p in (10, 50, 90)] for k, x in rel.items()}
    acc = np.ones(N, bool)
    for k in CAL_SET:
        acc &= np.abs(rel[k] - OBS[k]) <= TOL
    idx = np.where(acc)[0]
    post = {k: [float(np.percentile(rel[k][idx], p)) for p in (10, 50, 90)] for k in MIXES} if len(idx) else None
    hold_cov = float(np.mean(np.abs(rel[HOLDOUT][idx] - OBS[HOLDOUT]) <= TOL)) if len(idx) else None
    shifts = {}
    if len(idx):
        for p in ("kf", "nu", "bpack", "o_org", "eta_w", "kap_o", "kap_s", "kimp", "k0", "acap", "S0"):
            shifts[p] = dict(prior_median=float(np.median(d[p])), post_median=float(np.median(d[p][idx])), post_p10=float(np.percentile(d[p][idx], 10)), post_p90=float(np.percentile(d[p][idx], 90)))
    out["cases"][case] = dict(n=N, n_accepted=int(len(idx)), accept_rate=float(len(idx) / N), prior_pred=prior, post_pred=post, holdout_in_tol=hold_cov, shifts=shifts,
                              prior_pred_in_tol={k: float(np.mean(np.abs(rel[k] - OBS[k]) <= TOL)) for k in MIXES})
    print("case", case, "accepted", len(idx), "%.0fs" % (time.time() - T0), flush=True)
    if len(idx) < 200:
        continue
    dsub = {k: v[idx] for k, v in d.items()}
    res = {}
    opc = m.evaluate(case, dsub, dict(path="OPC"))
    res["opc"] = dict(L=[float(np.percentile(opc["L"], p)) for p in (10, 50, 90)], E=[float(np.percentile(opc["E_gate"], p)) for p in (10, 50, 90)], fc28=float(np.median(opc["fc28"])))
    for r in (0.10, 0.25, 0.50):
        for match in (False, True):
            sh = m.evaluate(case, dsub, dict(path="P1", r=r, d_um=30.0, match=match), ref=opc)
            res["r%d%s" % (round(100 * r), "m" if match else "u")] = m.summarize(sh, opc)
    # the same configurations on the prior draws (same size subsample) for a like-for-like comparison
    dp = {k: v[:6000] for k, v in d.items()}; dp["kf"] = np.zeros(6000); opcp = m.evaluate(case, dp, dict(path="OPC")); resp = {}
    for r in (0.10, 0.25, 0.50):
        for match in (False, True):
            resp["r%d%s" % (round(100 * r), "m" if match else "u")] = m.summarize(m.evaluate(case, dp, dict(path="P1", r=r, d_um=30.0, match=match), ref=opcp), opcp)
    out["cases"][case]["posterior_claims"] = res; out["cases"][case]["prior_claims"] = resp
    # carbon checks against the same paper (cement 0.91 kg/kg, shell powder 0.059 kg/kg measured lab value)
    if case == "A":
        ef = shell_ef(d, 55.0)
        out["shell_ef_model"] = dict(p10=float(np.percentile(ef, 10)), p50=float(np.median(ef)), p90=float(np.percentile(ef, 90)), published_lab=0.059)
        carbon = {}
        agg_w = 1390 * 0.00747 + 510 * 0.004 + 160 * 0.000344
        E0 = 350 * 0.91 + agg_w
        for k, mix in {"S15F15": (255, 95), "S30F15": (210, 140)}.items():
            Ek = mix[0] * 0.91 + mix[1] * 0.059 + agg_w
            carbon[k] = dict(published_E=round(Ek, 1), published_reduction_pct=round(100 * (Ek / E0 - 1), 1))
        carbon["control_E"] = round(E0, 1)
        # model prediction of carbon change at constant binder content with the published factors
        for k, r in (("S15F15", 95 / 350), ("S30F15", 140 / 350)):
            carbon[k]["model_reduction_pct_published_EF"] = round(100 * ((1 - r) * 350 * 0.91 + r * 350 * 0.059 + agg_w) / E0 - 100, 1)
        out["carbon_check"] = carbon
out["prior_ef_shell_note"] = "model shell EF at 55 um (kg CO2e/kg), prior draws"
json.dump(out, open("shell_out/calibration.json", "w"), indent=1)
print("done %.0fs" % (time.time() - T0))
