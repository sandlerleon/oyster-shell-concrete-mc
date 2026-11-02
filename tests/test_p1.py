"""Checks on the calibration, validation and primary results and on the original model (from the repository root: python -m pytest tests, or python tests/test_p1.py)."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
J = lambda *p: json.load(open(os.path.join(HERE, *p)))
CAL, KF0, VAL = J("shell_out", "calibration.json"), J("shell_out", "kf0_pred.json"), J("shell_out", "validation.json")
S1, BL, ALT = J("shell_cal", "stage1.json"), J("shell_cal", "blend.json"), J("shell_cal", "alternatives.json")


def test_original_model_is_rejected_by_the_scallop_data():
    for case in "ABC":
        assert KF0[case]["joint_cal_accept"] == 0
        assert all(v == 0.0 for v in KF0[case]["in_tol"].values())


def test_calibration_accepts_enough_draws_and_holds_out_mix():
    for case in "ABC":
        c = CAL["cases"][case]
        assert c["n_accepted"] > 5000 and c["holdout_in_tol"] > 0.95 and 0.3 < c["shifts"]["kf"]["post_median"] < 1.0


def test_study_level_kf_differs_in_sign_between_studies():
    est = VAL["kf_study_estimates"]
    assert est["KWON"] < 0 < est["ABBAS"] < est["ABIN"]
    assert VAL["kf_re"]["df"] == 3


def test_frozen_calibration_has_no_clear_advantage_on_independent_data():
    cmp_ = VAL["comparison"]
    a = cmp_["Abbas-calibrated posterior (frozen)"]["all"]["MAE"]; b = cmp_["Bolomey dilution only"]["all"]["MAE"]
    assert a > 0.5 * b and cmp_["Abbas-calibrated posterior (frozen)"]["KWON"]["bias"] > 0.3


def test_primary_results_use_the_pooled_distribution():
    assert S1["claims"]["r25"]["fc_rel"] > 0.8 and -23 < S1["claims"]["r25"]["dE_p50"] < -19


def test_joint_advantage_is_zero_against_blended_cement():
    for case in "ABC":
        for k in ("P1_r25u_vs_blend50", "P1_r25m_vs_blend50", "P1_r50m_vs_blend50"):
            assert BL["cases"][case][k]["p_adv_mat"] < 0.01


def test_porosity_coupling_raises_the_pavement_probability():
    assert ALT["claims"]["porosity_coupled"]["A"]["r25m"]["p_adv_mat"] > ALT["claims"]["strength_only"]["A"]["r25m"]["p_adv_mat"]


def test_leave_one_study_out_conclusions_are_stable():
    L = J("shell_cal", "loso_claims.json")
    pa = [v["claims"]["A"]["r25m"]["p_adv_mat"] for v in L.values()]; pc = [v["claims"]["C"]["r25m"]["p_adv_mat"] for v in L.values()]
    assert 0.1 < min(pa) and max(pa) < 0.3 and max(pc) < 0.06


def test_transport_scenarios_bracket_the_conclusion():
    t = ALT["transport"]
    assert t["chi0.0|coupled"]["A"]["r25m"]["p_adv_mat"] > 0.5 > t["chi1.0|strength_only"]["A"]["r25m"]["p_adv_mat"]


def test_baseline_recovery_of_original_model():
    import ShellConcreteMC_v2 as m
    v = m.verify()
    assert max(v["baseline_recovery_max_abs_diff"].values()) == 0.0 and v["carbon_closure_max_abs"] == 0.0


def test_kf_zero_reproduces_v2_strength_and_defaults_are_unchanged():
    import ShellConcreteMC_v2 as m2
    import ShellConcreteMC_v3 as m3
    d = m3.sample("A", 200, seed=1); B = np.full(200, 350.0); w = np.full(200, 160.0); T = np.full(200, 20.0)
    a = m2.fc28(d, B, w, 0.25, 30.0, T); d["kf"] = np.zeros(200)
    assert np.allclose(a, m3.fc28(d, B, w, 0.25, 30.0, T))
    assert np.allclose(m2.evaluate("A", m3.sample("A", 100, 2), dict(path="OPC"))["E_gate"], m3.evaluate("A", m3.sample("A", 100, 2), dict(path="OPC"))["E_gate"])


def test_blended_reference_cuts_carbon_and_extends_life():
    import ShellConcreteMC_v3 as m3
    d = m3.sample("A", 500, 3); o = m3.evaluate("A", d, dict(path="OPC")); b = m3.evaluate("A", d, dict(path="BLEND"))
    assert np.median(b["E_gate"] / o["E_gate"]) < 0.7 and np.median(b["L"]) > np.median(o["L"])


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print("ok", name)
