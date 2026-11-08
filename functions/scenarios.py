"""Run simulations for different clock configurations, error models, and
clock fractional-frequency uncertainties.
"""

import pandas as pd
from .constants import C2_G0, SIGMA_F
from .adjustment import error_cov, summarise


SHORT = {"Conventional": "conv", "Case 1": "c1", "Case 2": "c2", "Case 3": "c3"}
MODELS = {"random": "", "scale": "scale", "temp": "temp", "scale+temp": "scale+temp", "block": "block"}


def run_scenarios(A, sig_lev, sec, layouts, models=MODELS, sigma_f=SIGMA_F):
    """Evaluate all combinations of clock layouts, error models, and clock
    uncertainties.

    For each scenario, calculate the height uncertainty statistics using
    the selected estimator. Store the benchmark standard deviation profiles
    from the covariance-aware estimator for later plotting and comparison.

    Parameters
    ----------
    A : array-like
        Design matrix of the adjustment system.
    sig_lev : array-like
        Standard deviations of the levelling observations.
    sec : pandas.DataFrame
        Levelling section data used to construct the error covariance.
    layouts : dict
        Clock configurations for the conventional case and each clock case.
    models : dict
        Mapping between error-model labels and their corresponding model
        identifiers.
    sigma_f : sequence
        Clock fractional-frequency uncertainty values to evaluate.

    Returns
    -------
    pandas.DataFrame
        Summary statistics for each scenario, including RMSE, end-point
        uncertainty, worst benchmark uncertainty, and clock information.
    prof : dict
        Height standard deviation profiles from the covariance-aware
        estimator, indexed by case, error model, and clock uncertainty.
    """
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
    """Calculate benchmark height uncertainty profiles when systematic
    levelling errors are included.

    Evaluate each clock configuration using both the naive and
    covariance-aware estimators, and store the resulting standard
    deviation profiles for comparison.
    """
    out = {}
    for case in ["Case 1", "Case 2", "Case 3"]:
        for est in ["naive", "aware"]:
            Cx, _, _ = error_cov(A, sig_lev, model, layouts[case], sig_f, est, sec)
            out[(SHORT[case], est)] = summarise(Cx)[1]
    return out


def sensitivity_scale(A, sig_lev, sec, layouts, sig_kappas=(0, 5e-6, 10e-6, 20e-6, 30e-6), sig_f=1e-18):
    """Assess how rod-scale uncertainty affects the adjusted benchmark
    height uncertainties.

    Repeat the adjustment for each specified scale uncertainty and clock
    configuration, keeping the temperature-error model and clock
    uncertainty fixed. Return the summary statistics for each combination.
    """
    rows = []
    for sk in sig_kappas:
        for case, nodes in layouts.items():
            Cx, _, _ = error_cov(A, sig_lev, "scale+temp", nodes, sig_f, "aware", sec, sig_kappa=sk)
            rows.append({"sig_kappa_ppm": sk * 1e6, "case": case, **summarise(Cx)[0]})
    return pd.DataFrame(rows)


def sensitivity_delta(A, sec, layouts, delta_mm, sig_f=1e-18):
    """Evaluate the effect of levelling noise on the adjusted height
    uncertainties.

    Convert the specified section noise level from mm/√km to observation
    standard deviations using the section lengths. Run the adjustment
    for each clock configuration and return the resulting summary statistics.
    """
    import numpy as np
    sig = delta_mm * 1e-3 * np.sqrt(sec.L_km.values); rows = []
    for case, nodes in layouts.items():
        Cx, _, _ = error_cov(A, sig, "", nodes, sig_f, "aware", sec)
        rows.append({"case": case, "delta": delta_mm, **summarise(Cx)[0]})
    return pd.DataFrame(rows)