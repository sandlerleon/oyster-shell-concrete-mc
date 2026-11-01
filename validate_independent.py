# -*- coding: utf-8 -*-
"""Independent oyster-shell checks of the strength model with frozen parameters, comparison with simpler models, per-study filler efficiency
and a random-effects predictive distribution of k_f across studies.

    python validate_independent.py   ->  shell_out/validation.json
Data (28-day compressive strength relative to the control mix of the same study):
  ABBAS   Abbas and Kudukkan 2025 (scallop shell; calibration set of calibrate_abbas.py)
  KWON    Kwon and Wang 2024, Buildings 14:1578, Tables 1 and 3 (oyster shell powder, mortar, constant binder 568.16 kg/m3, w/b 0.50)
  ABIN    Abinaya et al., oyster shell, w/b 0.40, replacement 2.5-10 %: values as transcribed in Kim and Lee 2023, Buildings 13:3078, Table A1 (secondary source; binder content assumed)
  MEX     Mortar with untreated oyster-shell CaCO3, 10 wt% of cement mass, cement:sand 1:4, w/c 0.68 (Materials 19:1656, 2026, Table 2); sieve apertures as particle size
"""
import json
import time

import numpy as np

import os

import ShellConcreteMC_v3 as m

m.KF_SIZE_DEP = os.environ.get("KF_SIZE", "0") == "1"
OUTFILE = "shell_out/validation_size.json" if m.KF_SIZE_DEP else "shell_out/validation.json"

T0 = time.time()
TOL = 0.07
N = 300000


def mixes():
    S = {}
    S["ABBAS"] = dict(B=350.0, w=160.0, items=[("S0F15", 50 / 350, 90.0, 1.00), ("S15F15", 95 / 350, 64.0, (1 - 0.07) ** 2), ("S30F15", 140 / 350, 55.0, 0.75)], note="scallop shell, concrete")
    S["KWON"] = dict(B=568.16, w=284.08, items=[("OSP15", 0.15, 30.0, 32.1 / 49.2), ("OSP30", 0.30, 30.0, 21.1 / 49.2)], note="oyster shell powder, mortar; particle size not reported (30 um assumed)")
    S["ABIN"] = dict(B=400.0, w=160.0, items=[("r2.5", 0.025, 30.0, 40.0 / 38.0), ("r5", 0.05, 30.0, 42.0 / 38.0), ("r7.5", 0.075, 30.0, 40.0 / 38.0), ("r10", 0.10, 30.0, 39.5 / 38.0)],
                     note="oyster shell, concrete; binder content assumed 400 kg/m3, size 30 um assumed; secondary transcription")
    B, w = 2100.0 / 5.68, 0.68 * 2100.0 / 5.68
    S["MEX"] = dict(B=B, w=w, items=[("0.14 mm", 0.10, 140.0, 14.1 / 13.6), ("0.17 mm", 0.10, 170.0, 13.3 / 13.6), ("0.25 mm", 0.10, 250.0, 13.6 / 13.6), ("0.42 mm", 0.10, 420.0, 13.8 / 13.6), ("0.84 mm", 0.10, 840.0, 12.0 / 13.6)],
                    note="oyster-shell CaCO3, untreated, mortar 1:4, w/c 0.68; size = sieve aperture; 10 wt% interpreted as replacement")
    return S


def rel(d, B, w, r, d_um):
    n = len(d["B"]); Bv = np.full(n, B); wv = np.full(n, w); T = np.full(n, 20.0)
    return m.fc28(d, Bv, wv, r, d_um, T) / m.fc28(d, Bv, wv, 0.0, d_um, T)


def bolomey(r, k, wc_inv):
    """Dilution / k-value rule (Bolomey form): f ~ (c_eff/w - 0.5); relative to the control with c/w = wc_inv."""
    ceff = (1 - r) * wc_inv + k * r * wc_inv
    return (ceff - 0.5) / (wc_inv - 0.5)


S = mixes(); d = m.sample("A", N, seed=2024)
rng = np.random.default_rng(99); d["kf"] = rng.uniform(-1.5, 1.5, N)
out = {"tol": TOL, "datasets": {}, "per_study_kf": {}, "frozen": {}, "comparison": {}}
# ---- per-study ABC for k_f
R = {}
for s, v in S.items():
    R[s] = {nm: rel(d, v["B"], v["w"], r, du) for nm, r, du, _ in v["items"]}
    acc = np.ones(N, bool)
    for nm, r, du, obs in v["items"]:
        acc &= np.abs(R[s][nm] - obs) <= TOL
    idx = np.where(acc)[0]
    out["datasets"][s] = dict(note=v["note"], B=v["B"], w=v["w"], items=[dict(name=nm, r=r, d_um=du, observed=float(o)) for nm, r, du, o in v["items"]])
    if len(idx) >= 30:
        k = d["kf"][idx]
        out["per_study_kf"][s] = dict(n_accepted=int(len(idx)), median=float(np.median(k)), p10=float(np.percentile(k, 10)), p90=float(np.percentile(k, 90)), mean=float(k.mean()), sd=float(k.std(ddof=1)))
    else:
        # fall back to the best fitting k_f by least squares on a grid (no draw reaches the tolerance)
        grid = np.linspace(-1.5, 1.5, 61); sse = []
        dm = {kk: vv[:2000] for kk, vv in d.items()}
        for kf in grid:
            dm["kf"] = np.full(2000, kf); sse.append(sum(np.median((rel(dm, v["B"], v["w"], r, du) - o) ** 2) for _, r, du, o in v["items"]))
        out["per_study_kf"][s] = dict(n_accepted=int(len(idx)), best_grid_kf=float(grid[int(np.argmin(sse))]), min_sse=float(min(sse)))
    print(s, out["per_study_kf"][s], "%.0fs" % (time.time() - T0), flush=True)

# ---- frozen predictions of the Abbas-calibrated posterior and comparison models
ab = np.where(np.all([np.abs(R["ABBAS"][nm] - o) <= TOL for nm, _, _, o in S["ABBAS"]["items"] if nm != "S15F15"], axis=0))[0]
if len(ab) == 0:      # the formulation cannot reproduce the calibration set within the tolerance: report the per-study fits only
    out["note"] = "no draw reproduces the calibration mixes within the tolerance; frozen predictions not computed"
    json.dump(out, open(OUTFILE, "w"), indent=1); print("done (no accepted draws)"); raise SystemExit(0)
dpost = {kk: vv[ab] for kk, vv in d.items()}
dz = {kk: vv[:20000] for kk, vv in d.items()}; dz["kf"] = np.zeros(20000)
d25 = dict(dz); d25["kf"] = np.full(20000, 0.25)
rows = []
for s in ("KWON", "ABIN", "MEX"):
    v = S[s]
    for nm, r, du, o in v["items"]:
        post = rel(dpost, v["B"], v["w"], r, du); z = rel(dz, v["B"], v["w"], r, du); k25 = rel(d25, v["B"], v["w"], r, du)
        rows.append(dict(study=s, mix=nm, r=r, observed=float(o), post=[float(np.percentile(post, p)) for p in (5, 50, 95)], kf0=[float(np.percentile(z, p)) for p in (5, 50, 95)], kval25=[float(np.percentile(k25, p)) for p in (5, 50, 95)],
                         bolomey_dilution=float(bolomey(r, 0.0, v["B"] / v["w"])), bolomey_k25=float(bolomey(r, 0.25, v["B"] / v["w"]))))
out["frozen"]["rows"] = rows


def metrics(rows, key, pi=None):
    e = np.array([(r[key][1] if isinstance(r[key], list) else r[key]) - r["observed"] for r in rows])
    res = dict(n=len(e), MAE=float(np.abs(e).mean()), RMSE=float(np.sqrt((e ** 2).mean())), bias=float(e.mean()))
    if pi:
        res["coverage90"] = float(np.mean([(r[key][0] - TOL <= r["observed"] <= r[key][2] + TOL) for r in rows])); res["coverage90_no_tol"] = float(np.mean([(r[key][0] <= r["observed"] <= r[key][2]) for r in rows]))
    return res


for key, lab in (("post", "Abbas-calibrated posterior (frozen)"), ("kf0", "original inert-filler model (k_f = 0)"), ("kval25", "k-value rule in the same model (k_f = 0.25)"), ("bolomey_dilution", "Bolomey dilution only"), ("bolomey_k25", "Bolomey k-value 0.25")):
    for grp, sel in (("all", ("KWON", "ABIN", "MEX")), ("KWON", ("KWON",)), ("ABIN", ("ABIN",)), ("MEX", ("MEX",))):
        sub = [r for r in rows if r["study"] in sel]
        out["comparison"].setdefault(lab, {})[grp] = metrics(sub, key, pi=key in ("post", "kf0", "kval25"))

# ---- random-effects predictive distribution of k_f across studies (study-level posterior medians, reference prior)
ks = np.array([out["per_study_kf"][s].get("median", out["per_study_kf"][s].get("best_grid_kf")) for s in ("ABBAS", "KWON", "ABIN", "MEX")])
out["kf_study_estimates"] = dict(zip(("ABBAS", "KWON", "ABIN", "MEX"), map(float, ks)))


def re_pred(x):
    n = len(x); mu = x.mean(); sd = x.std(ddof=1)
    return dict(n=n, mu=float(mu), tau=float(sd), scale=float(sd * np.sqrt(1 + 1 / n)), df=n - 1)


out["kf_re"] = re_pred(ks)
loso = {}
for i, s in enumerate(("ABBAS", "KWON", "ABIN", "MEX")):
    keep = np.delete(ks, i); p = re_pred(keep)
    # predictive draws of k_f for the held-out study and resulting strength predictions for its mixes
    from scipy import stats
    kd = np.clip(p["mu"] + p["scale"] * stats.t.rvs(p["df"], size=20000, random_state=5 + i), -1.0, 1.5)
    dd = {kk: vv[:20000] for kk, vv in d.items()}; dd["kf"] = kd
    errs = []; cov = []
    for nm, r, du, o in S[s]["items"]:
        pr = rel(dd, S[s]["B"], S[s]["w"], r, du); q = [float(np.percentile(pr, x)) for x in (5, 50, 95)]
        errs.append(q[1] - o); cov.append(q[0] - TOL <= o <= q[2] + TOL)
    loso[s] = dict(kf_pred=p, MAE=float(np.mean(np.abs(errs))), bias=float(np.mean(errs)), coverage90=float(np.mean(cov)), n=len(errs))
out["loso"] = loso
json.dump(out, open(OUTFILE, "w"), indent=1)
print("done %.0fs" % (time.time() - T0))
