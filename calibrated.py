# -*- coding: utf-8 -*-
"""Calibrated variant of ShellConcreteMC_v3 for the full analysis.

`import calibrated as m` returns the v3 model module with the filler efficiency k_f added as an uncertain common input whose distribution is the random-effects
predictive distribution across the four published studies (shell_out/validation.json, written by validate_independent.py): a Student-t with location mu, scale
tau*sqrt(1 + 1/n) and n - 1 degrees of freedom, truncated to [-1, 1.5].   Environment variable KF_MODE: 're' (default), or 'zero' (k_f = 0, the original model).
"""
import json
import os
import sys

import numpy as np
from scipy import stats

import ShellConcreteMC_v3 as m

MODE = os.environ.get("KF_MODE", "re")
HERE = os.path.dirname(os.path.abspath(__file__))
V = json.load(open(os.path.join(HERE, "shell_out", "validation.json")))["kf_re"]
_old = m._ppf


def _ppf(spec, u):
    if spec[0] == "tre":
        return np.clip(spec[1] + spec[2] * stats.t.ppf(np.clip(u, 1e-9, 1 - 1e-9), spec[3]), spec[4], spec[5])
    return _old(spec, u)


m._ppf = _ppf
m.COMMON["kf"] = ("tre", V["mu"], V["scale"], V["df"], -1.0, 1.5) if MODE == "re" else ("uni", 0.0, 0.0)
m.LAB["kf"] = "Filler efficiency k_f"
sys.modules[__name__] = m
