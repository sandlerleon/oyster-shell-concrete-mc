#!/usr/bin/env python3
"""
ShellConcreteMC_v2 -- coupled hydration / transport / damage / carbon framework for
oyster-shell-derived cementitious materials in aggressive infrastructure environments.

Single importable file (numpy + scipy + matplotlib), same architecture as the author's MyUncle
framework: registries of uncertain inputs, closed-form solutions cross-checked numerically,
paired (common-random-number) comparison of each shell mix with its OPC reference, Latin
hypercube sampling, global (Sobol) sensitivity, falsification maps, one-call figure suite.

Pathways (NOT interchangeable)
  P1  uncalcined shell powder replacing a mass fraction r of cement (filler / diluent)
  P2  shell CaCO3 replacing a fraction rho of the limestone in clinker raw meal
  P3  shell-based alkali-activated binders: OUT OF SCOPE

All input ranges are literature-style placeholders, not measurements.
"""
import argparse, json, os
import numpy as np
from scipy import stats
from scipy.special import erfc
from scipy.stats import qmc

R_GAS = 8.314
KF_COUPLE_POROSITY = False     # v3 switch: strength-only efficiency (False) or efficiency that also refines the pores (True)
KF_SIZE_DEP = False             # v3 switch: efficiency independent of particle size (False) or scaled by fineness(d)/fineness(45 um) (True)
DT = 0.125
T_AX = np.arange(DT, 100.01, DT)
NT = len(T_AX)
H_CAP = 100.0
CO2_CALC = 44.01 / 100.09
RHO_C, RHO_F, RHO_W = 3150.0, 2710.0, 1000.0

COMMON = {
    "k0": ("logn", 0.26, 0.20), "nh": ("uni", 1.5, 2.5), "Ea_h": ("uni", 30e3, 45e3), "acap": ("uni", 0.90, 0.97),
    "S0": ("logn", 78.0, 0.10), "kap_s": ("uni", 1.8, 3.0), "nu": ("uni", 0.0, 1.5), "bpack": ("uni", 0.0, 0.30),
    "o_org": ("tri", 0.005, 0.02, 0.04), "eta_w": ("uni", 0.3, 0.9), "kap_o": ("logu", 5.0, 200.0), "kimp": ("uni", 0.0, 0.05),
    "lamD": ("logn", 1.0, 0.40), "chi_D": ("uni", 0.0, 1.0), "dzeta": ("uni", -0.3, 0.3), "m0": ("tri", 0.12, 0.20, 0.30), "omega": ("uni", 0.0, 0.30),
    "kapc": ("uni", 2.0, 20.0), "UD": ("uni", 35e3, 50e3), "Ccrit": ("tri", 0.6, 1.0, 1.6), "tprop": ("tri", 3.0, 6.0, 12.0),
    "eps0": ("tri", 400.0, 600.0, 900.0), "tau_sh": ("tri", 20.0, 35.0, 60.0), "KR": ("uni", 0.2, 0.6), "psi": ("uni", 0.4, 0.7),
    "d_cr": ("uni", 0.03, 0.10), "f_env": ("uni", 0.05, 0.25),
    "a_air": ("tri", 0.02, 0.05, 0.07), "S_crit": ("uni", 0.80, 0.92), "d1": ("logu", 3e-4, 3e-3), "q_ft": ("uni", 1.0, 3.0),
    "var_t": ("tri", 0.85, 1.0, 1.15), "n_traf": ("logu", 5e4, 1e6), "Dfail": ("tri", 0.2, 0.3, 0.4),
    "Qcem": ("uni", 350.0, 480.0), "psiT": ("uni", 0.5, 0.9), "Rres": ("uni", 0.3, 0.7), "alphaT": ("uni", 8e-6, 12e-6),
    "Csat": ("uni", 15.0, 22.0), "xi_L": ("uni", 0.3, 0.7), "kK": ("logu", 0.03, 1.0), "x_lim_L": ("tri", 30.0, 50.0, 100.0),
    "Dratio_Ca": ("uni", 0.5, 1.5), "xiN": ("uni", 0.4, 0.8), "dc": ("logu", 50.0, 300.0), "DH": ("uni", 3.0, 8.0),
    "x_lim_A": ("tri", 10.0, 20.0, 30.0), "k_ab": ("tri", 0.05, 0.30, 1.0), "b_ab": ("uni", 1.5, 2.5), "w_lim": ("tri", 15.0, 25.0, 40.0),
    "EFcem": ("tri", 0.75, 0.85, 0.95), "mCaCO3": ("tri", 1.0, 1.1, 1.2), "EFagg": ("tri", 0.004, 0.007, 0.012),
    "EFmisc": ("tri", 5.0, 10.0, 20.0), "e_dry": ("tri", 20.0, 60.0, 150.0), "E_ref": ("tri", 25.0, 45.0, 90.0),
    "E_cr": ("tri", 5.0, 10.0, 20.0), "EF_el": ("tri", 0.1, 0.4, 0.8), "EF_th": ("tri", 0.20, 0.23, 0.30),
    "dist_sh": ("tri", 30.0, 150.0, 600.0), "EF_tr": ("tri", 0.07, 0.10, 0.15), "c_land": ("tri", 0.0, 0.005, 0.02),
    "EF_lq": ("tri", 0.003, 0.006, 0.012), "kiln_pen": ("tri", 0.0, 0.005, 0.02), "p_sh": ("uni", 0.92, 0.98),
}
CASES = {
    "A": dict(name="Highway rigid pavement", B=("tri", 320.0, 350.0, 400.0), wb=("tri", 0.38, 0.42, 0.46), Cs=("tri", 2.0, 5.0, 9.0),
              cover=("tri", 50.0, 70.0, 100.0), Tenv=("tri", 2.0, 10.0, 18.0), RH=("uni", 55.0, 90.0), a_f=("uni", 0.1, 0.5),
              n_ft=("tri", 10.0, 40.0, 90.0), Tcure=("tri", 5.0, 20.0, 35.0), Econ=("tri", 15.0, 30.0, 55.0)),
    "B": dict(name="Dam (mass concrete)", B=("tri", 180.0, 230.0, 280.0), wb=("tri", 0.45, 0.50, 0.60), Cs=("tri", 0.0, 0.0001, 0.0002),
              cover=("tri", 50.0, 70.0, 100.0), Tenv=("tri", 4.0, 10.0, 16.0), RH=("uni", 90.0, 100.0), a_f=("uni", 0.2, 0.7),
              n_ft=("tri", 0.0, 15.0, 50.0), Tcure=("tri", 10.0, 15.0, 25.0), Econ=("tri", 10.0, 20.0, 40.0), grad=("tri", 5.0, 20.0, 60.0)),
    "C": dict(name="Marine culvert / riverbed", B=("tri", 340.0, 370.0, 420.0), wb=("tri", 0.38, 0.42, 0.46), Cs=("tri", 3.0, 6.0, 12.0),
              cover=("tri", 50.0, 65.0, 90.0), Tenv=("tri", 8.0, 16.0, 26.0), RH=("uni", 85.0, 100.0), a_f=("uni", 0.3, 0.9),
              n_ft=("tri", 0.0, 20.0, 60.0), Tcure=("tri", 5.0, 20.0, 35.0), Econ=("tri", 10.0, 25.0, 45.0)),
}
LOAD = {0: 0.28, 1: 0.38, 2: 0.46}
ACID = {0: 0.03, 1: 0.3, 2: 3.0}
CSMULT = {0: 0.5, 1: 1.0, 2: 2.0}
LEACHMULT = {0: 0.5, 1: 1.0, 2: 1.5}


def _ppf(spec, u):
    k = spec[0]; u = np.clip(u, 1e-9, 1 - 1e-9)
    if k == "tri":
        lo, mode, hi = spec[1:4]
        return stats.triang.ppf(u, (mode - lo) / (hi - lo), loc=lo, scale=hi - lo)
    if k == "uni":
        return spec[1] + (spec[2] - spec[1]) * u
    if k == "logn":
        return stats.lognorm.ppf(u, spec[2], scale=spec[1])
    if k == "logu":
        return np.exp(np.log(spec[1]) + (np.log(spec[2]) - np.log(spec[1])) * u)
    raise ValueError(k)


def _central(spec):
    k = spec[0]
    return {"tri": lambda: spec[2], "uni": lambda: 0.5 * (spec[1] + spec[2]), "logn": lambda: spec[1],
            "logu": lambda: float(np.sqrt(spec[1] * spec[2]))}[k]()


def specs_for(case):
    return {**COMMON, **{k: v for k, v in CASES[case].items() if k != "name"}}


def names_for(case):
    return list(specs_for(case))


def sample(case, n, seed=1, U=None, central=False):
    specs = specs_for(case); names = list(specs)
    if central:
        return {k: np.array([_central(s)]) for k, s in specs.items()}
    if U is None:
        U = qmc.LatinHypercube(d=len(names), seed=seed).random(n)
    return {nm: _ppf(specs[nm], U[:, i]) for i, nm in enumerate(names)}


def fineness(d_um):
    return float(np.clip((30.0 / d_um) ** 0.5, 0.0, 2.0))


def hydration(d, B, w, r, d_um, T_C, tdays, organic=True):
    """alpha(t) closed form for d(beta)/dt = -k beta^n, beta = alpha_max - alpha, alpha(0)=0."""
    c = np.maximum((1 - r) * B, 1e-6); wc = w / c
    amax = np.minimum(d["acap"], wc / 0.42)
    q = d["o_org"] * (1 - d["eta_w"]) * r / max(1 - r, 1e-3) if organic else 0.0
    A = np.exp(d["Ea_h"] / R_GAS * (1 / 293.15 - 1 / (T_C + 273.15)))
    k = d["k0"] * A * (1 + d["nu"] * fineness(d_um) * r) / (1 + d["kap_o"] * q)
    n = d["nh"][:, None]; t = np.atleast_1d(tdays)[None, :]
    base = amax[:, None] ** (1 - n) + (n - 1) * k[:, None] * t
    return amax[:, None] - base ** (1 / (1 - n)), amax


def mixture(B, w, r):
    return (1 - r) * B / RHO_C, r * B / RHO_F, w / RHO_W


def fc_series(d, alpha, vc, vw, r, d_um, rho_p2=0.0):
    kf = d.get("kf")                                     # v3: optional filler efficiency (k-value): the strength model counts c + kf*s as cementing mass
    if kf is not None and r > 0:
        if KF_SIZE_DEP:
            kf = kf * fineness(d_um) / fineness(45.0)
        vc = vc * np.maximum(1.0 + kf * r / max(1.0 - r, 1e-3), 0.0)      # effective cementing mass cannot be negative
    X = 2.06 * alpha * vc[:, None] / (vw[:, None] + alpha * vc[:, None])
    pack = 1 + d["bpack"] * fineness(d_um) * r * (1 - r)
    S0e = d["S0"] * 0.742 ** (2.0 - d["kap_s"])          # keeps nominal OPC strength (w/c 0.45) fixed for any exponent
    return S0e[:, None] * X ** d["kap_s"][:, None] * pack[:, None] * (1 - d["kimp"] * rho_p2)[:, None]


def fc28(d, B, w, r, d_um, T_C, rho_p2=0.0, organic=True):
    a, _ = hydration(d, B, w, r, d_um, T_C, 28.0, organic)
    vc, vf, vw = mixture(B, w, r)
    return fc_series(d, a, vc, vw, r, d_um, np.full(len(B), rho_p2))[:, 0]


def match_binder(d, B, w, r, d_um, T_C, target, organic=True):
    """Binder content B' (water fixed) with fc28 equal to the OPC target. Searches only the range in which the
    model's 28-day strength increases with binder (up to its peak); returns (B', matched_flag)."""
    sc = np.linspace(1.0, 2.5, 31)
    F = np.stack([fc28(d, B * s_, w, r, d_um, T_C, 0.0, organic) for s_ in sc], axis=1)
    ipk = np.argmax(F, axis=1); fpk = F[np.arange(len(B)), ipk]
    matched = fpk >= target
    lo = B.copy(); hi = B * sc[ipk]
    for _ in range(35):
        mid = 0.5 * (lo + hi); f = fc28(d, mid, w, r, d_um, T_C, 0.0, organic)
        lo = np.where(f < target, mid, lo); hi = np.where(f < target, hi, mid)
    return np.where(matched, 0.5 * (lo + hi), B * sc[ipk]), matched


def _calibrate():
    cen = {k: np.array([_central(s)]) for k, s in COMMON.items()}
    wcs = np.linspace(0.30, 0.60, 13); phis = []; Ds = []
    for wc in wcs:
        B = np.array([350.0]); w = wc * B
        a, _ = hydration(cen, B, w, 0.0, 30.0, 20.0, 28.0)
        vc, vf, vw = mixture(B, w, 0.0)
        phis.append(max(float(((vw - 1.06 * a[:, 0] * vc) / (vc + vf + vw))[0]), 1e-3)); Ds.append(10 ** (-12.06 + 2.40 * wc))
    zeta, lnD0 = np.polyfit(np.log(phis), np.log(Ds), 1)
    pred = np.exp(np.polyval([zeta, lnD0], np.log(phis)))
    return dict(zeta=float(zeta), lnD0=float(lnD0), wc=wcs.tolist(), phi=phis, D_life365=Ds, D_model=pred.tolist(),
                max_ratio=float(np.max(np.maximum(pred / Ds, Ds / pred))))


CAL = _calibrate()


def D28_of(d, phi_p):
    return d["lamD"] * np.exp(CAL["lnD0"]) * np.maximum(phi_p, 1e-3) ** (CAL["zeta"] + d["dzeta"])


def sigmoid(x):
    return 1 / (1 + np.exp(-np.clip(x, -40, 40)))


def evaluate(case, d, cfg, ref=None, keep=False):
    path = cfg.get("path", "OPC"); r = cfg.get("r", 0.0) if path == "P1" else 0.0
    fs = cfg.get("fs", 0.5) if path == "BLEND" else 0.0          # v3: blended-cement reference, slag fraction of the binder (screening treatment, see manuscript)
    rho = cfg.get("rho", 0.0) if path == "P2" else 0.0
    d_um = cfg.get("d_um", 30.0); sev = cfg.get("sev", 1); load = cfg.get("load", 1); bbio = cfg.get("beta_bio", 0.0)
    N = len(d["B"]); B0 = d["B"].copy(); w = d["wb"] * B0; Tc = d["Tcure"]
    if path == "P1" and cfg.get("match", False) and r > 0:
        B, matched = match_binder(d, B0, w, r, d_um, Tc, fc28(d, B0, w, 0.0, d_um, Tc))
    else:
        B = B0; matched = np.ones(len(B0), bool)
    c = (1 - r) * B; s_mass = r * B
    vc, vf, vw = mixture(B, w, r)
    rho_arr = np.full(N, rho)
    a6, amax = hydration(d, B, w, r, d_um, Tc, np.array([1.0, 3.0, 7.0, 28.0, 90.0, 365.0]))
    ag, _ = hydration(d, B, w, r, d_um, Tc, T_AX * 365.25)
    fcg = fc_series(d, ag, vc, vw, r, d_um, rho_arr)
    fc28v = fc_series(d, a6[:, 3:4], vc, vw, r, d_um, rho_arr)[:, 0]
    Eg = 4700.0 * np.sqrt(fcg); frg = 0.62 * np.sqrt(fcg)
    a28 = a6[:, 3]
    vce = vc
    if KF_COUPLE_POROSITY and d.get("kf") is not None and r > 0:
        kfe = d["kf"] * (fineness(d_um) / fineness(45.0) if KF_SIZE_DEP else 1.0); vce = vc * np.maximum(1.0 + kfe * r / max(1.0 - r, 1e-3), 0.0)
    phi_cap = np.maximum(vw - 1.06 * a28 * vce, 0.002); phi_p = phi_cap / (vc + vf + vw)
    out = dict(matched=matched, fc28=fc28v, alpha=a6, phi_cap=phi_cap, phi_p=phi_p, B=B, c=c, s=s_mass, fr28=0.62 * np.sqrt(fc28v))

    h = lambda x: 1 - (x / 100.0) ** 3
    eps_u = d["eps0"] * 1e-6 * (2.06 * a28 * vc / 0.183) * h(d["RH"]) / h(50.0)
    td = np.array([3.0, 7.0, 14.0, 28.0, 56.0, 90.0, 180.0, 365.0])
    a8, _ = hydration(d, B, w, r, d_um, Tc, td)
    f8 = fc_series(d, a8, vc, vw, r, d_um, rho_arr)
    eps8 = eps_u[:, None] * td[None, :] / (td[None, :] + d["tau_sh"][:, None])
    I_sh = np.max(d["KR"][:, None] * 4700 * np.sqrt(f8) * eps8 * (1 - d["psi"][:, None]) / np.maximum(0.33 * np.sqrt(f8), 1e-3), axis=1)
    a5, _ = hydration(d, B, w, r, d_um, Tc, 5.0)
    dTad = c * d["Qcem"] / 2400.0 * a5[:, 0]
    fcs = np.maximum(fc28v, 1e-3)
    I_th = d["Rres"] * 0.7 * 4700 * np.sqrt(fcs) * d["alphaT"] * d["psiT"] * dTad * (1 - d["psi"]) / (0.33 * np.sqrt(fcs))
    out.update(I_sh=I_sh, I_th=I_th, dTad=dTad)
    D0 = np.zeros(N)
    if case == "A": D0 = d["d_cr"] * sigmoid((I_sh - 1.0) / 0.15)
    if case == "B": D0 = d["d_cr"] * sigmoid((I_th - 1.0) / 0.15)

    D_phi = D28_of(d, phi_p)
    D_wc = d["lamD"] * 10 ** (-12.06 + 2.40 * np.clip(w / np.maximum(c, 1e-6), 0.25, 1.5))     # cement-only w/c mapping (pessimistic for filler)
    D28 = np.exp((1 - d["chi_D"]) * np.log(D_phi) + d["chi_D"] * np.log(D_wc))
    fT = np.exp(d["UD"] / R_GAS * (1 / 293.15 - 1 / (d["Tenv"] + 273.15)))
    m = np.minimum(d["m0"] + d["omega"] * r + 0.4 * fs / 0.7, 0.6)      # Life-365 treatment: slag raises the age exponent, D28 at equal w/cm unchanged
    Dt = (D28 * fT)[:, None] * (28.0 / (28.0 + T_AX * 365.25)[None, :]) ** m[:, None]
    ysec = 365.25 * 86400.0; tsec = T_AX * ysec
    Dch = np.zeros((N, NT)); Dab = np.zeros((N, NT)); x_acid = np.zeros((N, NT)); x_leach = np.zeros((N, NT))
    Df = d["Dfail"]
    if case == "B":
        Dca = Dt * d["Dratio_Ca"][:, None]
        NCa = np.maximum(d["xi_L"] * 11.6 * a28 * c, 1.0)
        Cs_ = d["Csat"] * LEACHMULT[sev]
        xd = np.sqrt(2 * np.cumsum(Dca * DT * ysec, axis=1) * Cs_[:, None] / NCa[:, None])
        xa = (d["kK"] * d["Dratio_Ca"] * D28 * fT * d["grad"] * Cs_ / NCa)[:, None] * tsec[None, :]
        x_leach = np.sqrt(xd ** 2 + xa ** 2) * 1000.0
        Dch = Df[:, None] * x_leach / d["x_lim_L"][:, None]
    if case == "C":
        chi = 1.0 / (1.0 + (d_um / d["dc"]) ** 2) if path == "P1" else np.zeros(N)
        Nacid = np.maximum(2.0 * d["xiN"] * 11.6 * a28 * c + 2.0 * chi * (s_mass * d["p_sh"]) / 0.10009, 1.0)
        x_acid = np.sqrt(2 * np.cumsum(Dt * d["DH"][:, None] * DT * ysec, axis=1) * ACID[sev] / Nacid[:, None]) * 1000.0
        Dch = Df[:, None] * x_acid / d["x_lim_A"][:, None]
        wear = np.cumsum(d["k_ab"][:, None] * (30.0 / np.maximum(fcg, 1.0)) ** d["b_ab"][:, None] * DT, axis=1)
        Dab = Df[:, None] * wear / d["w_lim"][:, None]

    fr_ref = ref["fr28"] if ref is not None else out["fr28"]
    sigma_traf = LOAD[load] * d["var_t"] * fr_ref if case == "A" else np.zeros(N)
    Sw = (phi_cap + d["a_f"] * d["a_air"]) / (phi_cap + d["a_air"])
    ft_rate = d["n_ft"] * d["d1"] * (np.maximum(Sw - d["S_crit"], 0.0) / (1 - d["S_crit"])) ** d["q_ft"]
    sig_env = np.zeros((N, NT))
    if case == "A":
        tg = T_AX * 365.25
        sig_env = d["f_env"][:, None] * d["KR"][:, None] * Eg * (eps_u[:, None] * tg[None, :] / (tg[None, :] + d["tau_sh"][:, None])) * (1 - d["psi"][:, None])
    Dm = D0.copy(); Dtot = np.zeros((N, NT)); Dmech_t = np.zeros((N, NT))
    for j in range(NT):
        D_now = Dm + Dch[:, j] + Dab[:, j]
        amp = 1.0 / np.clip(1.0 - D_now, 0.1, 1.0)               # Kachanov-Lemaitre effective-stress amplification
        rate = ft_rate.copy()
        if case == "A":
            S = np.clip(amp * (sigma_traf + sig_env[:, j]) / np.maximum(frg[:, j], 1e-3), 0.0, 1.6)
            rate = rate + d["n_traf"] / 10.0 ** np.maximum((1 - S) / (0.0685 * 0.9), -3.0)
        Dm = Dm + DT * rate
        Dmech_t[:, j] = Dm; Dtot[:, j] = Dm + Dch[:, j] + Dab[:, j]

    def first_cross(F, thr):
        reach = F >= thr[:, None]; j = np.argmax(reach, axis=1); anyr = reach.any(axis=1); jm = np.maximum(j - 1, 0)
        f1, f2 = F[np.arange(N), jm], F[np.arange(N), j]
        fr = np.clip((thr - f1) / np.maximum(f2 - f1, 1e-12), 0, 1)
        return np.where(anyr & (j > 0), T_AX[jm] + fr * DT, np.where(anyr, 0.1, np.inf))
    L_D = first_cross(Dtot, Df)
    L_cl = np.full(N, np.inf)
    if case in ("A", "C"):
        g = 1.0 + d["kapc"][:, None] * np.clip(Dmech_t, 0, 1)
        intD = np.cumsum(Dt * g * DT * ysec, axis=1)
        Cc = (d["Cs"] * CSMULT[sev])[:, None] * erfc((d["cover"] * 1e-3)[:, None] / (2 * np.sqrt(np.maximum(intD, 1e-30))))
        t_init = first_cross(Cc, d["Ccrit"]); L_cl = t_init + d["tprop"]; out["t_init"] = t_init
    L_raw = np.minimum(L_cl, L_D); L = np.minimum(L_raw, H_CAP)
    mode = np.where(L_cl <= L_D, 1, 2); mode = np.where(L_raw >= H_CAP, 0, mode)
    out.update(D100=Dtot[:, -1], Dfail=Df, L=L, L_cl=L_cl, L_D=L_D, mode=mode, D28=D28, m=m, Dtot=Dtot if keep else None,
               x_acid50=x_acid[:, int(round(50.0 / DT)) - 1], x_leach50=x_leach[:, int(round(50.0 / DT)) - 1])

    EF_sh = ((d["e_dry"] * d["EF_th"] + (d["E_ref"] * (45.0 / d_um) ** 0.5) * d["EF_el"]) / 1000.0 + d["dist_sh"] * d["EF_tr"] / 1000.0 - d["c_land"])
    EFcem = d["EFcem"].copy()
    if path == "BLEND":
        EFcem = (1.0 - fs) * EFcem + fs * cfg.get("EF_slag", 0.07)
    if path == "P2":
        sh_proc = (d["e_dry"] * d["EF_th"] + d["E_cr"] * d["EF_el"]) / 1000.0 + d["dist_sh"] * d["EF_tr"] / 1000.0 - d["c_land"] + d["kiln_pen"]
        calc = rho * d["mCaCO3"]
        EFcem = EFcem + calc * (sh_proc - d["EF_lq"]) - bbio * calc * CO2_CALC
    agg = 2400.0 - B - w
    E_bin = c * EFcem; E_sh = s_mass * EF_sh if path == "P1" else np.zeros(N)
    E_gate = E_bin + E_sh + agg * d["EFagg"] + d["EFmisc"]; E_tot = E_gate + d["Econ"]
    out.update(E_gate=E_gate, E_tot=E_tot, E_ann=E_tot / np.maximum(L, 1.0),
               parts=dict(binder=E_bin, shell=E_sh, aggregate=agg * d["EFagg"], water_admix=d["EFmisc"], construction=d["Econ"]))
    return out


def paired(case, d, cfg, ref=None):
    opc = ref if ref is not None else evaluate(case, d, dict(path="OPC", sev=cfg.get("sev", 1), load=cfg.get("load", 1)))
    return evaluate(case, d, cfg, ref=opc), opc


def noninferior(sh, opc, tol=0.0):
    """Service-life non-inferiority. Where the OPC reference survives the 100-yr horizon (censored), the shell mix must also survive
    and its damage index at the horizon must not exceed the reference by more than 5 % of the failure damage."""
    surv = opc["L"] >= H_CAP
    return np.where(surv, (sh["L"] >= H_CAP) & (sh["D100"] - opc["D100"] <= 0.05 * opc["Dfail"]), sh["L"] >= (1 - tol) * opc["L"])


def advantage(sh, opc, tol=0.0, annual=False, material=0.0):
    k = "E_ann" if annual else "E_gate"
    return (sh[k] < (1 - material) * opc[k]) & noninferior(sh, opc, tol)


def summarize(sh, opc):
    q = lambda a: [float(np.percentile(a, p)) for p in (10, 50, 90)]
    return dict(L=q(sh["L"]), L_opc=q(opc["L"]), dL=q(sh["L"] - opc["L"]), E=q(sh["E_gate"]), E_opc=q(opc["E_gate"]),
                dE_pct=q(100 * (sh["E_gate"] / opc["E_gate"] - 1)), Eann=q(sh["E_ann"]), fc28=q(sh["fc28"]), fc28_opc=q(opc["fc28"]),
                p_adv=float(np.mean(advantage(sh, opc))), p_adv_ann=float(np.mean(advantage(sh, opc, annual=True))),
                p_adv_tol5=float(np.mean(advantage(sh, opc, tol=0.05))), p_adv_mat=float(np.mean(advantage(sh, opc, material=0.05))),
                p_crack_sh=float(np.mean(sh["I_sh"] >= 1)), p_crack_th=float(np.mean(sh["I_th"] >= 1)), p_crack_th_opc=float(np.mean(opc["I_th"] >= 1)), p_crack_sh_opc=float(np.mean(opc["I_sh"] >= 1)), p_E=float(np.mean(sh["E_gate"] < opc["E_gate"])),
                p_L=float(np.mean(noninferior(sh, opc))), matched=float(np.mean(sh["matched"])), modes={str(k): float(np.mean(sh["mode"] == k)) for k in (0, 1, 2)})


# ---------------------------------------------------------------- verification
def _fd_check():
    from scipy.linalg import solve_banded
    D28, m, Cs = 6e-12, 0.25, 6.0; ysec = 365.25 * 86400; t_end = 40.0; nt = 4000; dt = t_end * ysec / nt
    nx = 600; x = np.linspace(0, 0.5, nx + 1); dx = x[1] - x[0]; Cc = np.zeros(nx + 1); Cc[0] = Cs
    for it in range(nt):
        D = D28 * (28.0 / (28.0 + (it + 0.5) * dt / ysec * 365.25)) ** m; lam = D * dt / dx ** 2
        ab = np.zeros((3, nx - 1)); ab[0, 1:] = -lam / 2; ab[1, :] = 1 + lam; ab[2, :-1] = -lam / 2
        rhs = Cc[1:-1] + lam / 2 * (Cc[2:] - 2 * Cc[1:-1] + Cc[:-2]); rhs[0] += lam / 2 * Cs
        Cc[1:-1] = solve_banded((1, 1), ab, rhs)
    ts = np.linspace(0, t_end, 40001)[1:]; Dtt = D28 * (28.0 / (28.0 + ts * 365.25)) ** m
    intD = np.trapezoid(Dtt, ts * ysec)
    err = float(np.max(np.abs(Cc - Cs * erfc(x / (2 * np.sqrt(intD))))) / Cs); assert err < 5e-3, err
    return dict(max_rel_err_vs_Cs=err)


def _dt_check():
    global DT, T_AX, NT
    d = sample("A", 400, 9); a = evaluate("A", d, dict(path="OPC")); old = (DT, T_AX, NT)
    try:
        DT = old[0] / 2; T_AX = np.arange(DT, 100.01, DT); NT = len(T_AX); b = evaluate("A", d, dict(path="OPC"))
    finally:
        DT, T_AX, NT = old
    rel = np.abs(b["L"] - a["L"]) / np.maximum(b["L"], 1.0)
    return dict(median_rel_diff=float(np.median(rel)), p95_rel_diff=float(np.percentile(rel, 95)))


def verify():
    res = {}; d = sample("A", 2000, 5); opc = evaluate("A", d, dict(path="OPC"))
    z1 = evaluate("A", d, dict(path="P1", r=0.0, d_um=25.0), ref=opc); z2 = evaluate("A", d, dict(path="P2", rho=0.0), ref=opc)
    res["baseline_recovery_max_abs_diff"] = dict(L_P1=float(np.max(np.abs(z1["L"] - opc["L"]))), E_P1=float(np.max(np.abs(z1["E_gate"] - opc["E_gate"]))),
                                                 L_P2=float(np.max(np.abs(z2["L"] - opc["L"]))), E_P2=float(np.max(np.abs(z2["E_gate"] - opc["E_gate"]))))
    assert max(res["baseline_recovery_max_abs_diff"].values()) < 1e-9
    dd = {k: v[:6] for k, v in d.items()}; B = dd["B"]; w = dd["wb"] * B; td = np.array([1.0, 7.0, 28.0, 90.0, 365.0])
    a, _ = hydration(dd, B, w, 0.2, 25.0, dd["Tcure"], td)
    amx = np.minimum(dd["acap"], w / (0.8 * B) / 0.42); q = dd["o_org"] * (1 - dd["eta_w"]) * 0.2 / 0.8
    k = dd["k0"] * np.exp(dd["Ea_h"] / R_GAS * (1 / 293.15 - 1 / (dd["Tcure"] + 273.15))) * (1 + dd["nu"] * fineness(25.0) * 0.2) / (1 + dd["kap_o"] * q)
    err = 0.0
    for i in range(6):
        al = 0.0; steps = 20000; h = 365.0 / steps; tt = []; vals = []
        f = lambda x: k[i] * max(amx[i] - x, 0) ** dd["nh"][i]
        for s_ in range(steps):
            k1 = f(al); k2 = f(al + .5 * h * k1); k3 = f(al + .5 * h * k2); k4 = f(al + h * k3); al += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            tt.append((s_ + 1) * h); vals.append(al)
        err = max(err, float(np.max(np.abs(np.interp(td, tt, vals) - a[i]))))
    res["hydration_closed_form_vs_rk4_max_abs_err"] = err; assert err < 1e-3, err
    res["chloride_erfc_vs_finite_difference"] = _fd_check(); res["damage_step_convergence"] = _dt_check()
    e = evaluate("A", d, dict(path="P1", r=0.25, d_um=25.0), ref=opc)
    clos = float(np.max(np.abs(sum(e["parts"].values()) - e["E_tot"]))); assert clos < 1e-9; res["carbon_closure_max_abs"] = clos
    res["calcination_stoichiometry_kgCO2_per_kgCaCO3"] = CO2_CALC
    res["life365_calibration"] = dict(zeta=CAL["zeta"], max_ratio_model_to_life365=CAL["max_ratio"])
    tgt = fc28(d, d["B"], d["wb"] * d["B"], 0.0, 25.0, d["Tcure"]); Bm, mt = match_binder(d, d["B"], d["wb"] * d["B"], 0.25, 25.0, d["Tcure"], tgt)
    ferr = np.abs(fc28(d, Bm, d["wb"] * d["B"], 0.25, 25.0, d["Tcure"]) / tgt - 1)
    added = mt & (Bm > d["B"] * (1 + 1e-9))      # v3: with k_f > 0 the shell mix can already exceed the target at the original binder content (no binder is added or removed)
    res["strength_matching_max_rel_err_matched_draws"] = float(np.max(ferr[added])) if added.any() else 0.0; res["strength_matching_feasible_fraction_r25"] = float(np.mean(mt)); assert res["strength_matching_max_rel_err_matched_draws"] < 5e-3
    res["status"] = "all assertions passed"
    return res


def qualitative_checks(n=3000, seed=3):
    d = sample("A", n, seed); base = evaluate("A", d, dict(path="OPC")); rs = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]
    fcs, phis = [], []
    for r in rs:
        e = evaluate("A", d, dict(path="P1", r=r, d_um=25.0), ref=base); fcs.append(e["fc28"]); phis.append(e["phi_cap"])
    fcs, phis = np.array(fcs), np.array(phis); res = {}
    res["strength_monotone_decrease_prob"] = float(np.mean(np.all(np.diff(fcs, axis=0) < 0, axis=0)))
    res["capillary_porosity_increase_prob_0_to_50"] = float(np.mean(phis[-1] > phis[0]))
    res["relative_strength_28d_median"] = {f"{int(r*100)}%": float(np.median(fcs[i] / fcs[0])) for i, r in enumerate(rs)}
    gel = []
    for r in (0.0, 0.20):
        e = evaluate("A", d, dict(path="P1", r=r, d_um=25.0), ref=base); gel.append(2.06 * e["alpha"][:, 3] * (1 - r) * e["B"] / RHO_C)
    res["shrinkage_proxy_decrease_prob_0_to_20"] = float(np.mean(gel[1] < gel[0]))
    d2 = {k: v.copy() for k, v in d.items()}; d2["eta_w"] = np.zeros(n); d2["nu"] = np.zeros(n)
    e0 = evaluate("A", d2, dict(path="OPC")); e1 = evaluate("A", d2, dict(path="P1", r=0.2, d_um=25.0), ref=e0)
    res["early_age_retardation_prob_unwashed_no_filler_effect"] = float(np.mean(e1["alpha"][:, 0] < e0["alpha"][:, 0]))
    res["early_age_acceleration_prob_default"] = float(np.mean(evaluate("A", d, dict(path="P1", r=0.2, d_um=25.0), ref=base)["alpha"][:, 0] > base["alpha"][:, 0]))
    return res


def sobol(case, cfg, n=1024, seed=11, chunk=6000, n_boot=100):
    names = names_for(case); dim = len(names); rng = np.random.default_rng(seed)
    A = qmc.Sobol(d=2 * dim, scramble=True, seed=seed).random(n); UA, UB = A[:, :dim], A[:, dim:]

    def run(U):
        L, E = [], []
        for s in range(0, len(U), chunk):
            sh, opc = paired(case, sample(case, 0, U=U[s:s + chunk]), cfg)
            L.append(sh["L"] - opc["L"]); E.append(100 * (sh["E_gate"] / opc["E_gate"] - 1))
        return np.concatenate(L), np.concatenate(E)
    yA, yB = run(UA), run(UB); store = {}
    for i in range(dim):
        U = UA.copy(); U[:, i] = UB[:, i]; store[i] = run(U)
    out = {}
    for oi, key in enumerate(("dL", "dE")):
        fA, fB = yA[oi], yB[oi]; V = np.var(np.concatenate([fA, fB])) + 1e-30; st, s1, stc = [], [], []
        for i in range(dim):
            fAB = store[i][oi]; st.append(float(np.mean((fA - fAB) ** 2) / (2 * V))); s1.append(float(1 - np.mean((fB - fAB) ** 2) / (2 * V)))
            bs = []
            for _ in range(n_boot):
                ix = rng.integers(0, n, n); v = np.var(np.concatenate([fA[ix], fB[ix]])) + 1e-30; bs.append(np.mean((fA[ix] - fAB[ix]) ** 2) / (2 * v))
            stc.append(np.percentile(bs, [2.5, 97.5]).tolist())
        out[key] = dict(names=names, ST=st, S1=s1, ST_ci=stc, sum_S1=float(np.sum(s1)))
    return out


R_GRID = [0.0, 0.10, 0.25, 0.50, 0.75, 1.0]; D_GRID = [10.0, 30.0, 100.0, 200.0]; RHO_GRID = [0.25, 0.50, 0.75, 1.0]; BIO = [0.0, 0.5, 1.0]


def study(case, n=3000, seed=21, sev=1, load=1):
    d = sample(case, n, seed); opc = evaluate(case, d, dict(path="OPC", sev=sev, load=load)); rows = []
    for match in (False, True):
        for r in R_GRID[1:]:
            for dum in D_GRID:
                sh = evaluate(case, d, dict(path="P1", r=r, d_um=dum, match=match, sev=sev, load=load), ref=opc)
                rows.append(dict(path="P1", r=r, d_um=dum, match=match, **summarize(sh, opc)))
    for rho in RHO_GRID:
        for b in BIO:
            sh = evaluate(case, d, dict(path="P2", rho=rho, beta_bio=b, sev=sev, load=load), ref=opc)
            rows.append(dict(path="P2", rho=rho, beta_bio=b, **summarize(sh, opc)))
    q = lambda a: [float(np.percentile(a, p)) for p in (10, 50, 90)]
    return dict(case=case, rows=rows, opc=dict(L=q(opc["L"]), E=q(opc["E_gate"]), fc28=float(np.median(opc["fc28"])), modes={str(k): float(np.mean(opc["mode"] == k)) for k in (0, 1, 2)}))


def severity_load_scan(case, cfg, n=3000, seed=31):
    d = sample(case, n, seed); out = {}
    combos = {"A": [(1, 0), (1, 1), (1, 2), (0, 1), (2, 1)], "B": [(s, 1) for s in range(3)], "C": [(s, 1) for s in range(3)]}[case]
    for sev, load in combos:
        opc = evaluate(case, d, dict(path="OPC", sev=sev, load=load)); out[f"sev{sev}_load{load}"] = summarize(evaluate(case, d, dict(cfg, sev=sev, load=load), ref=opc), opc)
    return out


# ---------------------------------------------------------------- figures
C_OPC, C_P1, C_P2, C_GR = "#555555", "#1b6ca8", "#c8553d", "#888888"
LAB = dict(k0="Hydration rate", nh="Reaction order", Ea_h="Hydration Ea", acap="Max. hydration", S0="Gel strength S0", kap_s="Strength exponent", nu="Filler nucleation", bpack="Packing bonus",
           o_org="Shell organics", eta_w="Washing efficiency", kap_o="Organic retardation", kimp="Impurity penalty", lamD="D28 factor", chi_D="D28 model form", dzeta="Porosity exponent",
           m0="Age exponent m", omega="m shift", kapc="Crack\u2013diffusion coupling", UD="Diffusion Ea", Ccrit="Critical chloride", tprop="Propagation period",
           eps0="Shrinkage strain", tau_sh="Shrinkage time", KR="Restraint", psi="Creep relaxation", d_cr="Early crack damage", f_env="Curling factor", a_air="Air content",
           S_crit="Critical saturation", d1="FT damage/cycle", q_ft="FT exponent", var_t="Load variability", n_traf="Traffic cycles/yr",
           Dfail="Failure damage", Qcem="Heat of hydration", psiT="Thermal drop", Rres="Thermal restraint", alphaT="Thermal expansion", Csat="Leaching saturation",
           xi_L="Leachable Ca", kK="Permeability factor", x_lim_L="Leaching limit", Dratio_Ca="Ca diffusivity", xiN="Acid-accessible Ca", dc="Calcite access size",
           DH="H+ diffusivity", x_lim_A="Acid limit", k_ab="Abrasion rate", b_ab="Abrasion exponent", w_lim="Wear limit", EFcem="Cement EF", mCaCO3="Raw-meal CaCO3",
           EFagg="Aggregate EF", EFmisc="Water/admixture EF", e_dry="Shell drying energy", E_ref="Grinding energy", E_cr="Crushing energy", EF_el="Electricity EF",
           EF_th="Heat EF", dist_sh="Shell haul distance", EF_tr="Truck EF", c_land="Landfill credit", EF_lq="Limestone EF", kiln_pen="Kiln penalty", p_sh="Shell purity",
           B="Binder content", wb="w/b", Cs="Surface chloride", cover="Cover", Tenv="Service temperature", RH="Relative humidity", a_f="Air-void filling",
           n_ft="Freeze\u2013thaw cycles/yr", Tcure="Curing temperature", Econ="Construction emissions", grad="Hydraulic gradient")


def _style():
    import matplotlib as mpl
    mpl.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25, "savefig.dpi": 220, "savefig.bbox": "tight"})


def fig_architecture(path):
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    fig, ax = plt.subplots(figsize=(10, 4.4)); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 44)
    def box(t, x, y, w, h, fc, ec):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3", fc=fc, ec=ec, lw=1.2)); ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=7.5)
    def arr(x0, y0, x1, y1): ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", color="#444", lw=1.2))
    box("P1  uncalcined shell powder\nreplaces cement fraction r\n(filler / diluent)", 1, 31, 25, 11, "#e8f0f8", C_P1)
    box("P2  shell CaCO$_3$ replaces limestone\nin clinker raw meal (fraction $\\rho$)", 1, 17, 25, 11, "#f4e8e4", C_P2)
    box("P3  shell-based alkali-activated\nbinders: OUT OF SCOPE", 1, 3, 25, 11, "#eeeeee", C_GR)
    box("Hydration\nporosity, strength, stiffness,\nearly-age cracking indices", 32, 29, 28, 13, "#e8f0f8", C_P1)
    box("Transport\nchloride (Fick, $D(t)$);\nacid and leaching fronts", 32, 14, 28, 12, "#e8f0f8", C_P1)
    box("Damage\nfatigue, freeze\u2013thaw, thermal,\nshrinkage, abrasion, chemical", 64, 29, 28, 13, "#e8f0f8", C_P1)
    box("Service life $L$\n(corrosion, damage limit)", 64, 14, 28, 12, "#e8f0f8", C_P1)
    box("Carbon accounting $E$: binder + shell + aggregate + construction $-$ credits\n(cradle-to-gate per m$^3$; annualised by $L$)", 32, 1, 60, 10, "#eaf5ea", "#2e7d32")
    arr(26.5, 36, 31.5, 36); arr(26.5, 22, 31.5, 22); arr(46, 28.6, 46, 26.4); arr(60.5, 36, 63.5, 36); arr(60.5, 20, 63.5, 20)
    arr(78, 28.6, 78, 26.4); arr(46, 13.6, 46, 11.4); arr(78, 13.6, 78, 11.4); arr(26.5, 8, 31.5, 6)
    fig.savefig(path); plt.close(fig)


def fig_hydration(path):
    import matplotlib.pyplot as plt
    d = sample("A", 1500, 7); ta = np.logspace(0, np.log10(365), 60); B = np.full(1500, 350.0); w = 0.42 * B
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.7))
    for r, col in ((0.0, C_OPC), (0.10, "#2e86c1"), (0.25, C_P1), (0.50, C_P2)):
        a, _ = hydration(d, B, w, r, 30.0, np.full(1500, 20.0), ta)
        axs[0].plot(ta, np.median(a, 0), color=col, label=f"r = {int(r*100)}%"); axs[0].fill_between(ta, np.percentile(a, 10, 0), np.percentile(a, 90, 0), color=col, alpha=0.12)
    axs[0].set_title("(a) Replacement fraction (30 \u00b5m, 20 \u00b0C)", fontsize=9); axs[0].set_ylabel("Hydration degree $\\alpha$")
    for T, col in ((5, "#2e86c1"), (20, C_OPC), (40, "#f39c12"), (60, C_P2)):
        a, _ = hydration(d, B, w, 0.25, 30.0, np.full(1500, float(T)), ta); axs[1].plot(ta, np.median(a, 0), color=col, label=f"{T} \u00b0C")
    axs[1].set_title("(b) Temperature (r = 25%)", fontsize=9)
    for dum, col in ((10, "#2e86c1"), (30, C_P1), (100, "#f39c12"), (200, C_P2)):
        a, _ = hydration(d, B, w, 0.25, float(dum), np.full(1500, 20.0), ta); axs[2].plot(ta, np.median(a, 0), color=col, label=f"{dum} \u00b5m")
    axs[2].set_title("(c) Particle size (r = 25%)", fontsize=9)
    for ax in axs: ax.set_xscale("log"); ax.set_xlabel("Curing time (d)"); ax.legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def fig_transport(path, n=3000):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.7)); rs = [0.0, 0.10, 0.25, 0.50]
    for ax, case, key, lab in ((axs[0], "A", "D28", "$D_{28}$ (m$^2$/s)"), (axs[1], "C", "x_acid50", "Acid front at 50 yr (mm)"), (axs[2], "B", "x_leach50", "Leaching front at 50 yr (mm)")):
        d = sample(case, n, 5); sev = 2 if case == "C" else 1; opc = evaluate(case, d, dict(path="OPC", sev=sev))
        for match, col, ls in ((False, C_P1, "-"), (True, C_P2, "--")):
            med, lo, hi = [], [], []
            for r in rs:
                e = opc if r == 0 else evaluate(case, d, dict(path="P1", r=r, d_um=30.0, match=match, sev=sev), ref=opc); v = e[key]
                med.append(np.median(v)); lo.append(np.percentile(v, 10)); hi.append(np.percentile(v, 90))
            ax.plot([100 * r for r in rs], med, ls, color=col, marker="o", label="same w/b" if not match else "strength-matched"); ax.fill_between([100 * r for r in rs], lo, hi, color=col, alpha=0.10)
        ax.set_xlabel("Shell replacement of cement (%)"); ax.set_ylabel(lab)
    axs[0].set_yscale("log"); axs[0].set_title("(a) Chloride diffusivity (Case A)", fontsize=9); axs[1].set_title("(b) Acid, severe (Case C)", fontsize=9); axs[2].set_title("(c) Soft-water leaching (Case B)", fontsize=9)
    axs[0].legend(frameon=False, fontsize=7); fig.tight_layout(); fig.savefig(path); plt.close(fig)


def fig_life(path, n=3000):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.8), sharey=True)
    ecdf = lambda v: (np.sort(v), np.arange(1, len(v) + 1) / len(v))
    for ax, case in zip(axs, "ABC"):
        d = sample(case, n, 9); opc = evaluate(case, d, dict(path="OPC")); x, y = ecdf(opc["L"]); ax.step(x, y, color=C_OPC, lw=2, label="OPC")
        for r, col in ((0.10, "#2e86c1"), (0.25, C_P1), (0.50, C_P2)):
            x, y = ecdf(evaluate(case, d, dict(path="P1", r=r, d_um=30.0), ref=opc)["L"]); ax.step(x, y, color=col, label=f"P1 r={int(r*100)}%")
        x, y = ecdf(evaluate(case, d, dict(path="P1", r=0.25, d_um=30.0, match=True), ref=opc)["L"]); ax.step(x, y, "--", color=C_P1, label="P1 r=25%, matched")
        ax.set_xlabel("Service life (yr, censored at 100)"); ax.set_title(CASES[case]["name"], fontsize=9)
    axs[0].set_ylabel("Cumulative probability"); axs[0].legend(frameon=False, fontsize=7); fig.tight_layout(); fig.savefig(path); plt.close(fig)


def pareto_mask(E, L):
    return np.array([not np.any((E <= E[i]) & (L >= L[i]) & ((E < E[i]) | (L > L[i]))) for i in range(len(E))])


def fig_pareto(path, studies):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.9))
    for ax, case in zip(axs, "ABC"):
        st = studies[case]; o = st["opc"]; pts = [(o["E"][1], o["L"][1], C_OPC, "o")]
        for r in st["rows"]:
            if r["path"] == "P1" and abs(r["d_um"] - 30.0) < 1e-9 and r["r"] <= 0.5: pts.append((r["E"][1], r["L"][1], "#2e86c1" if r["match"] else C_P1, "^" if r["match"] else "s"))
            if r["path"] == "P2" and r["beta_bio"] == 0.0 and r["rho"] == 1.0: pts.append((r["E"][1], r["L"][1], C_P2, "D"))
        E = np.array([p[0] / max(p[1], 1.0) for p in pts]); L = np.array([p[1] for p in pts]); keep = pareto_mask(E, L)
        for (e, l, col, mk), k, ea in zip(pts, keep, E): ax.scatter(ea, l, c=col if k else "none", edgecolors=col, marker=mk, s=40, zorder=3)
        pf = sorted([(ea, p[1]) for p, k, ea in zip(pts, keep, E) if k]); ax.step([p[0] for p in pf], [p[1] for p in pf], where="post", color="k", lw=0.8, alpha=0.6)
        for lab, col, mk in (("OPC", C_OPC, "o"), ("P1 same w/b", C_P1, "s"), ("P1 matched", "#2e86c1", "^"), ("P2 (no credit)", C_P2, "D")): ax.scatter([], [], c=col, marker=mk, label=lab)
        ax.set_xlabel("Embodied carbon per year of life (kg CO$_2$e/m$^3$/yr)"); ax.set_title(CASES[case]["name"], fontsize=9)
    axs[0].set_ylabel("Median service life (yr)"); axs[0].legend(frameon=False, fontsize=7); fig.tight_layout(); fig.savefig(path); plt.close(fig)


def fig_falsification(path, studies):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 3, figsize=(12.5, 6.4))
    for j, case in enumerate("ABC"):
        for i, match in enumerate((False, True)):
            ax = axs[i, j]; Z = np.zeros((len(D_GRID), len(R_GRID) - 1)); NM = np.zeros_like(Z, dtype=bool)
            for r in studies[case]["rows"]:
                if r["path"] == "P1" and r["match"] == match:
                    a_, b_ = D_GRID.index(r["d_um"]), R_GRID[1:].index(r["r"]); Z[a_, b_] = r["p_adv_mat"]; NM[a_, b_] = match and r["matched"] < 0.5
            im = ax.imshow(np.ma.masked_where(NM, Z), origin="lower", aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
            ax.imshow(np.where(NM, 1.0, np.nan), origin="lower", aspect="auto", cmap="Greys", vmin=0, vmax=3)
            ax.set_xticks(range(len(R_GRID) - 1), [f"{int(100*r)}" for r in R_GRID[1:]]); ax.set_yticks(range(len(D_GRID)), [int(x) for x in D_GRID])
            for a in range(Z.shape[0]):
                for b in range(Z.shape[1]): ax.text(b, a, "strength\nnot matchable" if NM[a, b] else f"{Z[a,b]:.2f}", ha="center", va="center", fontsize=5.5 if NM[a, b] else 7)
            ax.grid(False); ax.set_title(f"{CASES[case]['name']}\n{'strength-matched' if match else 'same w/b'}", fontsize=8)
            if i == 1: ax.set_xlabel("Replacement r (%)")
            if j == 0: ax.set_ylabel("Median particle size (\u00b5m)")
    fig.colorbar(im, ax=axs, label="P(\u22655 % carbon cut and non-inferior life)", shrink=0.8); fig.savefig(path); plt.close(fig)


def fig_p2(path, studies):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.6), sharey=True)
    for ax, case in zip(axs, "ABC"):
        for b, col in zip(BIO, ("#c8553d", "#f39c12", "#2e7d32")):
            ax.plot([100 * x for x in RHO_GRID], [r["dE_pct"][1] for r in studies[case]["rows"] if r["path"] == "P2" and r["beta_bio"] == b], "-o", color=col, label=f"biogenic credit {int(100*b)}%")
        ax.axhline(0, color="k", lw=0.8); ax.set_xlabel("Limestone replaced by shell in raw meal (%)"); ax.set_title(CASES[case]["name"], fontsize=9)
    axs[0].set_ylabel("Median change in embodied carbon (%)"); axs[0].legend(frameon=False, fontsize=7); fig.tight_layout(); fig.savefig(path); plt.close(fig)


def fig_sobol(path, sob):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 3, figsize=(12.5, 6.2))
    for j, case in enumerate("ABC"):
        for i, key in enumerate(("dL", "dE")):
            s = sob[case][key]; st = np.array(s["ST"]); ix = np.argsort(st)[-8:]; ax = axs[i, j]
            ax.barh([LAB.get(s["names"][k], s["names"][k]) for k in ix], st[ix], color=C_P1 if key == "dL" else C_P2)
            err = np.abs(np.array([[st[k] - s["ST_ci"][k][0], s["ST_ci"][k][1] - st[k]] for k in ix]).T)
            ax.errorbar(st[ix], range(len(ix)), xerr=err, fmt="none", ecolor="k", lw=0.8, capsize=2)
            ax.set_title(f"{CASES[case]['name']}: " + ("\u0394 service life" if key == "dL" else "\u0394 embodied carbon"), fontsize=8); ax.grid(axis="y", alpha=0)
            if i == 1: ax.set_xlabel("Total-order Sobol index")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def fig_verification(path, qc):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.6)); ax = axs[0]
    ax.loglog(CAL["D_life365"], CAL["D_model"], "o", color=C_P1); lim = [min(CAL["D_life365"]) * 0.7, max(CAL["D_life365"]) * 1.4]; ax.plot(lim, lim, "k--", lw=0.8)
    ax.set_xlabel("Life-365 $D_{28}$ (m$^2$/s)"); ax.set_ylabel("Model $D_{28}$ (m$^2$/s)"); ax.set_title("(a) OPC baseline, $w/c$ = 0.30\u20130.60", fontsize=9)
    d = sample("A", 1500, 3); B = np.full(1500, 350.0); wcs = np.linspace(0.30, 0.62, 12)
    f = np.array([np.percentile(fc28(d, B, wc * B, 0.0, 30.0, np.full(1500, 20.0)), [10, 50, 90]) for wc in wcs])
    axs[1].plot(wcs, f[:, 1], color=C_P1); axs[1].fill_between(wcs, f[:, 0], f[:, 2], alpha=0.15, color=C_P1)
    axs[1].set_xlabel("w/c"); axs[1].set_ylabel("28-day OPC strength (MPa)"); axs[1].set_title("(b) OPC strength vs w/c (P10\u2013P90)", fontsize=9)
    rs = [0, 5, 10, 15, 20, 30, 40, 50]; axs[2].plot(rs, [qc["relative_strength_28d_median"][f"{r}%"] for r in rs], "-o", color=C_P2)
    axs[2].set_xlabel("Shell replacement (%)"); axs[2].set_ylabel("Relative 28-day strength (median)"); axs[2].set_title("(c) Directional check vs published trend", fontsize=9)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def fig_carbon(path, n=3000):
    import matplotlib.pyplot as plt
    d = sample("A", n, 5); opc = evaluate("A", d, dict(path="OPC"))
    cf = [("OPC", opc)] + [(lab, evaluate("A", d, cfg, ref=opc)) for lab, cfg in (
        ("P1 10%", dict(path="P1", r=0.10, d_um=30.0)), ("P1 25%", dict(path="P1", r=0.25, d_um=30.0)), ("P1 25%\nmatched", dict(path="P1", r=0.25, d_um=30.0, match=True)),
        ("P1 50%", dict(path="P1", r=0.50, d_um=30.0)), ("P1 100%", dict(path="P1", r=1.0, d_um=30.0)),
        ("P2 100%\n0% credit", dict(path="P2", rho=1.0, beta_bio=0.0)), ("P2 100%\n100% credit", dict(path="P2", rho=1.0, beta_bio=1.0)))]
    keys = ["binder", "shell", "aggregate", "water_admix", "construction"]; cols = ["#555555", "#1b6ca8", "#c9b37e", "#9bb7d4", "#8c6d46"]
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 3.9)); bottom = np.zeros(len(cf))
    for k, col in zip(keys, cols):
        v = np.array([np.mean(c[1]["parts"][k]) for c in cf]); axs[0].bar(range(len(cf)), v, bottom=bottom, color=col, label=k.replace("_", "/")); bottom += v
    axs[0].set_xticks(range(len(cf)), [c[0] for c in cf], fontsize=6.5); axs[0].set_ylabel("kg CO$_2$e / m$^3$ (mean)"); axs[0].legend(frameon=False, fontsize=7, ncol=2)
    rs = np.linspace(0, 1, 11); med, lo, hi = [], [], []
    for r in rs:
        e = opc if r == 0 else evaluate("A", d, dict(path="P1", r=r, d_um=30.0), ref=opc); v = 100 * (e["E_gate"] / opc["E_gate"] - 1)
        med.append(np.median(v)); lo.append(np.percentile(v, 10)); hi.append(np.percentile(v, 90))
    axs[1].plot(100 * rs, med, "-o", color=C_P1); axs[1].fill_between(100 * rs, lo, hi, color=C_P1, alpha=0.15); axs[1].axhline(-80, color="#c8553d", ls="--", lw=1)
    axs[1].text(2, -77, "claimed 80% reduction", color="#c8553d", fontsize=7); axs[1].set_xlabel("Shell replacement of cement (%)"); axs[1].set_ylabel("Embodied-carbon change vs OPC (%)")
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def reproduce(outdir="out", n=3000, sobol_n=1024):
    os.makedirs(outdir, exist_ok=True); _style()
    ver = verify(); qc = qualitative_checks(); studies = {c: study(c, n) for c in "ABC"}
    base = dict(path="P1", r=0.25, d_um=30.0, match=False)
    scans = {c: severity_load_scan(c, base) for c in "ABC"}; scans_m = {c: severity_load_scan(c, dict(base, match=True)) for c in "ABC"}
    sob = {c: sobol(c, base, n=sobol_n) for c in "ABC"}
    d = sample("A", n, 5); opc = evaluate("A", d, dict(path="OPC")); claims = {}
    for r in (0.25, 0.5, 1.0):
        e = evaluate("A", d, dict(path="P1", r=r, d_um=30.0), ref=opc); em = evaluate("A", d, dict(path="P1", r=r, d_um=30.0, match=True), ref=opc)
        v = 100 * (e["E_gate"] / opc["E_gate"] - 1); vm = 100 * (em["E_gate"] / opc["E_gate"] - 1)
        claims[f"r{int(r*100)}"] = dict(dE_p50=float(np.median(v)), dE_p10=float(np.percentile(v, 10)), dE_p90=float(np.percentile(v, 90)), dE_matched_p50=float(np.median(vm)),
                                        fc_rel=float(np.median(e["fc28"] / opc["fc28"])), L_p50=float(np.median(e["L"])), L_opc_p50=float(np.median(opc["L"])))
    res = dict(n=n, sobol_n=sobol_n, verification=ver, qualitative=qc, calibration=CAL, studies=studies, scans=scans, scans_matched=scans_m,
               sobol={c: sob[c] for c in "ABC"}, claims=claims)
    fig_architecture(f"{outdir}/fig1_architecture.png"); fig_hydration(f"{outdir}/fig2_hydration.png"); fig_transport(f"{outdir}/fig3_transport.png")
    fig_life(f"{outdir}/fig4_life.png"); fig_pareto(f"{outdir}/fig5_pareto.png", studies); fig_falsification(f"{outdir}/fig6_falsification.png", studies)
    fig_p2(f"{outdir}/fig7_p2.png", studies); fig_sobol(f"{outdir}/fig8_sobol.png", sob); fig_verification(f"{outdir}/fig9_verification.png", qc); fig_carbon(f"{outdir}/fig10_carbon.png")
    json.dump(res, open(f"{outdir}/results.json", "w"), indent=1); return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verify", action="store_true"); ap.add_argument("--reproduce", action="store_true")
    ap.add_argument("--n", type=int, default=3000); ap.add_argument("--sobol-n", type=int, default=1024); ap.add_argument("--outdir", default="out")
    a = ap.parse_args()
    if a.verify: print(json.dumps(verify(), indent=1)); print(json.dumps(qualitative_checks(), indent=1)); print(json.dumps({k: CAL[k] for k in ("zeta", "lnD0", "max_ratio")}))
    elif a.reproduce: reproduce(a.outdir, a.n, a.sobol_n); print("done")
    else: ap.print_help()
