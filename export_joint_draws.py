# -*- coding: utf-8 -*-
"""Export joint draws of cement displacement, relative service life and carbon terms from the calibrated Paper 1 model for use in the allocation paper.

    python export_joint_draws.py <out.json>
For each case (A pavement, B dam, C culvert) and shell replacement r in (0.10, 0.25, 0.50), 30 um shell, strength-matched design, paired with its OPC reference on the same draw:
  delta  = cement avoided per tonne of shell used = (c_OPC - c_shell) / s_shell        (t/t; may be <= 0 when matching adds more cement than the shell replaces)
  rhoL   = L_shell / L_OPC                                                              (service-life ratio, both censored at 100 yr)
  E_sh, E_opc = cradle-to-gate carbon per m3 of shell concrete and of OPC concrete (kg CO2e/m3, shell processing and haul included, avoided-landfill credit included)
  s_kg   = shell mass per m3 of shell concrete (kg)
The joint draws preserve the correlations between displacement, life and carbon that arise from the shared mix design and exposure.
"""
import json
import sys

import numpy as np

import calibrated as m

OUT = sys.argv[1]
rows = {k: [] for k in ("case", "r", "delta", "rhoL", "E_sh", "E_opc", "s_kg", "kf", "L_sh", "L_opc", "c_sh", "c_opc")}
N = 3000
for ci, case in enumerate("ABC"):
    for r in (0.10, 0.25, 0.50):
        d = m.sample(case, N, seed=500 + 10 * ci + int(100 * r))
        opc = m.evaluate(case, d, dict(path="OPC"))
        sh = m.evaluate(case, d, dict(path="P1", r=r, d_um=30.0, match=True), ref=opc)
        ok = sh["matched"] & (sh["s"] > 1.0)
        delta = (opc["c"] - sh["c"]) / sh["s"]
        for k, v in dict(case=np.full(N, ci), r=np.full(N, r), delta=delta, rhoL=sh["L"] / opc["L"], E_sh=sh["E_gate"], E_opc=opc["E_gate"], s_kg=sh["s"], kf=d["kf"], L_sh=sh["L"], L_opc=opc["L"], c_sh=sh["c"], c_opc=opc["c"]).items():
            rows[k].extend([float(x) for x in np.asarray(v)[ok]])
rows["n"] = len(rows["case"])
json.dump(rows, open(OUT, "w"))
q = lambda a: [float(np.percentile(a, p)) for p in (10, 50, 90)]
print(rows["n"], "delta", q(rows["delta"]), "rhoL", q(rows["rhoL"]))
