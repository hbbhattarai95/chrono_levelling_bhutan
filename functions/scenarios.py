"""Scenario runs: clock layouts x error models x clock uncertainties."""
import pandas as pd
from .constants import C2_G0, SIGMA_F
from .adjustment import error_cov, summarise

SHORT = {"Conventional": "conv", "Case 1": "c1", "Case 2": "c2", "Case 3": "c3"}
MODELS = {"random": "", "scale": "scale", "temp": "temp", "scale+temp": "scale+temp", "block": "block"}


def run_scenarios(A, sig_lev, sec, layouts, models=MODELS, sigma_f=SIGMA_F):
    """Table of RMSE / end / worst for every layout, model, sigma_f and estimator.
    Also returns s.d. profiles (aware estimator) keyed by (case_short, model, sigma_f)."""
    rows, prof = [], {}
    for case, nodes in layouts.items():
        for mname, mod in models.items():
            for sf in (sigma_f if nodes is not None else [sigma_f[0]]):
                for est in (["aware", "naive"] if (mod and nodes is not None) else ["aware"]):
                    Cx, _, _ = error_cov(A, sig_lev, mod, nodes, sf, est, sec)
                    s, p = summarise(Cx)
                    rows.append({"case": case, "model": mname, "estimator": est,
                                 "n_clocks": 0 if nodes is None else len(nodes), "sig_f": sf,
                                 "sig_h_mm": C2_G0 * sf * 1e3, **s})
                    if est == "aware":
                        prof[(SHORT[case], mod, sf)] = p
    return pd.DataFrame(rows), prof


def systematic_profiles(A, sig_lev, sec, layouts, model="scale+temp", sig_f=1e-18):
    out = {}
    for case in ["Case 1", "Case 2", "Case 3"]:
        for est in ["naive", "aware"]:
            Cx, _, _ = error_cov(A, sig_lev, model, layouts[case], sig_f, est, sec)
            out[(SHORT[case], est)] = summarise(Cx)[1]
    return out


def sensitivity_scale(A, sig_lev, sec, layouts, sig_kappas=(0, 5e-6, 10e-6, 20e-6, 30e-6), sig_f=1e-18):
    rows = []
    for sk in sig_kappas:
        for case, nodes in layouts.items():
            Cx, _, _ = error_cov(A, sig_lev, "scale+temp", nodes, sig_f, "aware", sec, sig_kappa=sk)
            rows.append({"sig_kappa_ppm": sk * 1e6, "case": case, **summarise(Cx)[0]})
    return pd.DataFrame(rows)


def sensitivity_delta(A, sec, layouts, delta_mm, sig_f=1e-18):
    import numpy as np
    sig = delta_mm * 1e-3 * np.sqrt(sec.L_km.values); rows = []
    for case, nodes in layouts.items():
        Cx, _, _ = error_cov(A, sig, "", nodes, sig_f, "aware", sec)
        rows.append({"case": case, "delta": delta_mm, **summarise(Cx)[0]})
    return pd.DataFrame(rows)
