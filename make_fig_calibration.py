# -*- coding: utf-8 -*-
"""Figures for the out-of-sample test and the posterior-conditioned claims (reads shell_out/calibration.json).

    python make_fig_calibration.py <outdir>      ->  fig11_calibration.png, fig12_posterior.png
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else "shell_out"
C = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "shell_out", "calibration.json")))
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
NAVY, BLUE, GREY, RED = "#1f4e79", "#1f77b4", "#9aa5b1", "#c0392b"

# ---- Fig 11: predictive distributions against observations
cs = C["cases"]["A"]; tol = C["tol"]
mixes = [("S0F15", "S0F15\n(50 of 350 kg/m³ filler,\ncalibration)"), ("S15F15", "S15F15\n(95 kg/m³ shell,\nheld out)"), ("S30F15", "S30F15\n(140 kg/m³ shell,\ncalibration)")]
fig, ax = plt.subplots(1, 2, figsize=(10.2, 3.9), gridspec_kw={"width_ratios": [1.5, 1]})
a = ax[0]
for i, (k, lab) in enumerate(mixes):
    pr = cs["prior_pred"][k]; po = cs["post_pred"][k]; ob = C["obs"][k]
    a.add_patch(plt.Rectangle((i - 0.30, pr[0]), 0.26, pr[2] - pr[0], fc=GREY, ec="none", alpha=0.85, label="prior predictive P10–P90" if i == 0 else None))
    a.plot([i - 0.30, i - 0.04], [pr[1]] * 2, color="#445", lw=1.6)
    a.add_patch(plt.Rectangle((i + 0.04, po[0]), 0.26, max(po[2] - po[0], 0.004), fc=BLUE, ec="none", label="posterior predictive P10–P90" if i == 0 else None))
    a.plot([i + 0.04, i + 0.30], [po[1]] * 2, color="white", lw=1.4)
    a.add_patch(plt.Rectangle((i - 0.34, ob - tol), 0.68, 2 * tol, fc="none", ec=RED, ls=":", lw=1.0, label="observed ± tolerance" if i == 0 else None))
    a.plot([i], [ob], marker="D", color=RED, ms=7, ls="none", label="observed relative strength" if i == 0 else None)
a.set_xticks(range(3)); a.set_xticklabels([m[1] for m in mixes], fontsize=8); a.set_xlim(-0.6, 2.6); a.set_ylim(0.4, 1.15)
a.set_ylabel("28-day strength relative to control"); a.grid(axis="y", alpha=0.3); a.legend(loc="lower left", fontsize=7.5, frameon=False)
a.set_title("(a) Strength of scallop-shell concrete, constant binder", fontsize=9.5, loc="left", color=NAVY)
b = ax[1]
sh = cs["shifts"]; names = [("kf", "filler efficiency kf"), ("nu", "nucleation ν"), ("bpack", "packing β_p"), ("kap_o", "retardation κ_o"), ("eta_w", "washing η_w"), ("kap_s", "strength exponent κ")]
prior_rng = {"kf": (0.0, 1.2)}
for j, (k, lab) in enumerate(names):
    s = sh[k]; lo, hi, md = s["post_p10"], s["post_p90"], s["post_median"]; pm = s["prior_median"]
    b.plot([lo / pm, hi / pm], [j, j], color=BLUE, lw=5, solid_capstyle="butt"); b.plot([md / pm], [j], marker="o", color="white", mec=NAVY, ms=5, zorder=3)
b.axvline(1.0, color="#445", lw=1, ls="--"); b.set_yticks(range(len(names))); b.set_yticklabels([n[1] for n in names], fontsize=8); b.invert_yaxis(); b.set_xscale("log"); b.set_xticks([0.1, 0.3, 1, 3, 10]); b.set_xticklabels(["0.1", "0.3", "1", "3", "10"]); b.minorticks_off()
b.set_xlabel("posterior P10–P90 / prior median"); b.grid(axis="x", alpha=0.3)
b.set_title("(b) What the data constrain", fontsize=9.5, loc="left", color=NAVY)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig11_calibration.png"), dpi=200, facecolor="white"); plt.close(fig)

# ---- Fig 12: prior-model versus posterior-conditioned P(joint advantage)
keys = [("r10u", "10%\nsame w/b"), ("r25u", "25%\nsame w/b"), ("r25m", "25%\nmatched"), ("r50u", "50%\nsame w/b"), ("r50m", "50%\nmatched")]
titles = {"A": "(a) Pavement", "B": "(b) Dam", "C": "(c) Marine culvert"}
fig, ax = plt.subplots(1, 3, figsize=(10.4, 3.5), sharey=True)
for a, case in zip(ax, "ABC"):
    cc = C["cases"][case]
    pri = [cc["prior_claims"][k]["p_adv_mat"] for k, _ in keys]; pos = [cc["posterior_claims"][k]["p_adv_mat"] for k, _ in keys]
    x = range(len(keys))
    a.bar([i - 0.2 for i in x], pri, 0.38, color=GREY, label="original model (kf = 0)")
    a.bar([i + 0.2 for i in x], pos, 0.38, color=BLUE, label="calibrated (posterior)")
    for i, (u, v) in enumerate(zip(pri, pos)):
        a.text(i - 0.2, u + 0.012, "%.2f" % u, ha="center", fontsize=6.5, color="#445"); a.text(i + 0.2, v + 0.012, "%.2f" % v, ha="center", fontsize=6.5, color=NAVY)
    a.set_xticks(list(x)); a.set_xticklabels([l for _, l in keys], fontsize=7.5); a.set_title(titles[case], fontsize=9.5, loc="left", color=NAVY); a.grid(axis="y", alpha=0.3); a.set_ylim(0, 0.7)
ax[0].set_ylabel("P(≥5% carbon cut and non-inferior life)"); ax[0].legend(fontsize=7.5, frameon=False, loc="upper right")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig12_posterior.png"), dpi=200, facecolor="white"); plt.close(fig)
print("ok")
