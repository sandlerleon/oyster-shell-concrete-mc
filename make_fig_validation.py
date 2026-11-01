# -*- coding: utf-8 -*-
"""Figures: independent validation / per-study filler efficiency, and the ablation of the claims across the three model variants.
    python make_fig_validation.py <outdir>  ->  fig12_validation.png, fig13_ablation.png"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else "shell_cal"
V = json.load(open("shell_out/validation.json")); CAL = json.load(open("shell_out/calibration.json")); S = json.load(open("shell_cal/stage1.json"))
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
NAVY, BLUE, TEAL, ORANGE, GREY, RED = "#1f4e79", "#1f77b4", "#2a9d8f", "#e08a1e", "#9aa5b1", "#c0392b"

fig, ax = plt.subplots(1, 2, figsize=(10.4, 4.0), gridspec_kw={"width_ratios": [1.1, 1]})
a = ax[0]; rows = V["frozen"]["rows"]; col = {"KWON": RED, "ABIN": ORANGE, "MEX": TEAL}
for r in rows:
    o = r["observed"]; p = r["post"]
    a.plot([o, o], [p[0], p[2]], color=col[r["study"]], lw=1.2, alpha=0.8); a.plot([o], [p[1]], "o", color=col[r["study"]], ms=5)
    a.plot([o], [r["kf0"][1]], "x", color=GREY, ms=5)
lim = (0.3, 1.2); a.plot(lim, lim, "k--", lw=0.9); a.set_xlim(*lim); a.set_ylim(*lim); a.set_aspect("equal")
a.plot([], [], "o", color=NAVY, label="Abbas-calibrated, frozen (median, 5–95%)"); a.plot([], [], "x", color=GREY, label="original model (k_f = 0)")
for k, l in (("KWON", "Kwon & Wang 2024"), ("ABIN", "Abinaya et al. (secondary)"), ("MEX", "Mortar, Mexico (2026)")): a.plot([], [], "s", color=col[k], label=l)
a.set_xlabel("observed relative 28-day strength"); a.set_ylabel("predicted relative strength"); a.legend(fontsize=6.8, frameon=False, loc="upper left"); a.grid(alpha=0.3)
a.set_title("(a) Independent oyster-shell data", fontsize=9.5, loc="left", color=NAVY)
b = ax[1]; est = V["per_study_kf"]; names = [("ABBAS", "Abbas & Kudukkan\n(scallop, concrete)"), ("KWON", "Kwon & Wang\n(oyster, mortar)"), ("ABIN", "Abinaya et al.\n(oyster, secondary)"), ("MEX", "Mortar 2026\n(oyster CaCO₃)")]
for i, (k, lab) in enumerate(names):
    e = est[k]
    if "median" in e:
        b.plot([e["p10"], e["p90"]], [i, i], color=BLUE, lw=6, solid_capstyle="butt", alpha=0.85); b.plot([e["median"]], [i], "o", mfc="white", mec="k", ms=6, zorder=3)
re = V["kf_re"]; from scipy import stats
lo, hi = [re["mu"] + re["scale"] * stats.t.ppf(p, re["df"]) for p in (0.1, 0.9)]
b.plot([max(lo, -1.0), min(hi, 1.5)], [4.3, 4.3], color=ORANGE, lw=6, solid_capstyle="butt"); b.plot([re["mu"]], [4.3], "o", mfc="white", mec="k", ms=6, zorder=3)
b.axvline(0, color="k", lw=0.8); b.set_yticks(range(5)); b.set_yticklabels([n[1] for n in names] + ["random-effects\npredictive"], fontsize=8); b.invert_yaxis(); b.set_xlabel("effective strength-efficiency coefficient k_f (median and P10–P90)"); b.grid(axis="x", alpha=0.3)
b.set_title("(b) k_f differs between studies", fontsize=9.5, loc="left", color=NAVY)
fig.tight_layout(); fig.savefig(OUT + "/fig12_validation.png", dpi=200, facecolor="white"); plt.close(fig)

# ablation of the joint-advantage probability across variants
cfgs = [("A", "r25u", "pavement\n25% same w/b"), ("A", "r25m", "pavement\n25% matched"), ("B", "r25m", "dam\n25% matched"), ("C", "r25u", "culvert\n25% same w/b"), ("C", "r25m", "culvert\n25% matched")]
def pooled(case, key):
    r = float(key[1:3]) / 100; mt = key.endswith("m")
    return [x for x in S["studies"][case]["rows"] if x["path"] == "P1" and abs(x["r"] - r) < 1e-9 and x["d_um"] == 30 and x["match"] == mt][0]["p_adv_mat"]
fig, a = plt.subplots(figsize=(8.8, 3.5)); x = np.arange(len(cfgs)); w = 0.26
v0 = [CAL["cases"][c]["prior_claims"][k]["p_adv_mat"] for c, k, _ in cfgs]; v1 = [CAL["cases"][c]["posterior_claims"][k]["p_adv_mat"] for c, k, _ in cfgs]; v2 = [pooled(c, k) for c, k, _ in cfgs]
for off, v, c, lab in ((-w, v0, GREY, "original inert filler (k_f = 0)"), (0, v1, BLUE, "calibrated on Abbas & Kudukkan only"), (w, v2, TEAL, "pooled across four studies (primary)")):
    a.bar(x + off, v, w, color=c, label=lab)
    for xi, vi in zip(x + off, v): a.text(xi, vi + 0.01, "%.2f" % vi, ha="center", fontsize=6.8)
a.set_xticks(x); a.set_xticklabels([c[2] for c in cfgs], fontsize=8); a.set_ylabel("P(≥5% carbon cut and non-inferior life)"); a.legend(fontsize=7.5, frameon=False); a.grid(axis="y", alpha=0.3); a.set_ylim(0, 0.45)
fig.tight_layout(); fig.savefig(OUT + "/fig13_ablation.png", dpi=200, facecolor="white"); plt.close(fig)
print("ok")
