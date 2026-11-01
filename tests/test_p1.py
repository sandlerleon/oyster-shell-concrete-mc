"""Checks on the calibration results and the original model (run from the repository root: python -m pytest tests, or python tests/test_p1.py)."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAL = json.load(open(os.path.join(HERE, "shell_out", "calibration.json")))
KF0 = json.load(open(os.path.join(HERE, "shell_out", "kf0_pred.json")))


def test_original_model_is_rejected_by_the_data():
    for case in "ABC":
        assert KF0[case]["joint_cal_accept"] == 0
        assert all(v == 0.0 for v in KF0[case]["in_tol"].values())


def test_calibration_accepts_enough_draws_and_holds_out_mix():
    for case in "ABC":
        c = CAL["cases"][case]
        assert c["n_accepted"] > 5000
        assert c["holdout_in_tol"] > 0.95
        assert 0.3 < c["shifts"]["kf"]["post_median"] < 1.0


def test_posterior_predictive_within_tolerance_of_observations():
    for m, obs in CAL["obs"].items():
        lo, med, hi = CAL["cases"]["A"]["post_pred"][m]
        assert abs(med - obs) <= CAL["tol"] + 1e-9, (m, med, obs)


def test_baseline_recovery_of_original_model():
    import ShellConcreteMC_v2 as m
    v = m.verify()
    assert max(v["baseline_recovery_max_abs_diff"].values()) == 0.0
    assert v["carbon_closure_max_abs"] == 0.0


def test_kf_zero_reproduces_v2_strength():
    import numpy as np
    import ShellConcreteMC_v2 as m2
    import ShellConcreteMC_v3 as m3
    d = m3.sample("A", 200, seed=1); B = np.full(200, 350.0); w = np.full(200, 160.0); T = np.full(200, 20.0)
    a = m2.fc28(d, B, w, 0.25, 30.0, T)
    d["kf"] = np.zeros(200)
    b = m3.fc28(d, B, w, 0.25, 30.0, T)
    assert np.allclose(a, b)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print("ok", name)
