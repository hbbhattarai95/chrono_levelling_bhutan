"""Monte Carlo simulations to verify the exact propagation of observation
uncertainties through the adjustment.
"""

import numpy as np
import pandas as pd
from .adjustment import error_cov, summarise


def mc_verify(A, sig_lev, sec, nodes, sig_f, model, n_mc, rng, d_mm):
    """Compare formal uncertainty estimates with Monte Carlo results.

    Calculate the formal RMSE from the adjusted covariance matrix, then
    estimate the RMSE using two simulation approaches:

    1. Gaussian Monte Carlo: generate random errors from the full observation
       covariance matrix.
    2. Bootstrap Monte Carlo: resample the random levelling component from
       the observed normalised misclosures, while retaining Gaussian
       simulation for the remaining error sources.

    The bootstrap resampling uses the observed misclosures normalised by
    twice the square root of their section lengths. Resampled values are
    scaled to the corresponding section lengths and assigned random signs.

    Returns
    -------
    formal : float
        Formal RMSE from the adjusted covariance matrix, in millimetres.
    gauss : float
        RMSE estimated using Gaussian Monte Carlo simulations, in millimetres.
    boot : float
        RMSE estimated using bootstrap Monte Carlo simulations, in millimetres.
    """
    Cx, G, S = error_cov(A, sig_lev, model, nodes, sig_f, "aware", sec)
    formal = summarise(Cx)[0]["rmse_mm"]

    # Generate Gaussian observation errors using the full covariance matrix.
    eg = rng.standard_normal((n_mc, len(S))) @ np.linalg.cholesky(S).T
    gauss = 1e3 * np.sqrt(np.mean((eg @ G.T) ** 2))

    # Separate the non-random covariance contribution from the levelling
    # observation covariance before generating the bootstrap simulations.
    m = A.shape[0]
    S_rest = S.copy(); S_rest[:m, :m] -= np.diag(sig_lev ** 2)
    w, V = np.linalg.eigh(S_rest); w = np.clip(w, 0, None)
    er = rng.standard_normal((n_mc, len(S))) @ (V * np.sqrt(w)).T

    # Resample the observed normalised levelling misclosures, apply random
    # signs, and scale them according to the corresponding section lengths.
    zk = d_mm * 1e-3 / (2 * np.sqrt(sec.L_km.values))
    er[:, :m] += rng.choice(zk, (n_mc, m)) * rng.choice([-1.0, 1.0], (n_mc, m)) * np.sqrt(sec.L_km.values)
    boot = 1e3 * np.sqrt(np.mean((er @ G.T) ** 2))

    return formal, gauss, boot


def run_verification(A, sig_lev, sec, layouts, d_mm, n_mc=5000, rng=None):
    """Verify the formal adjustment uncertainties for all clock scenarios.

    Run Gaussian and bootstrap Monte Carlo simulations for the conventional
    levelling configuration and each clock layout. Both random-only and
    systematic levelling error models are considered.

    The results are collected in a DataFrame containing the formal RMSE,
    the two Monte Carlo RMSE estimates, and their ratios relative to the
    formal RMSE.
    """
    rng = rng or np.random.default_rng(20261104)
    rows = []
    for case in ["Conventional", "Case 1", "Case 2", "Case 3"]:
        for mod in ["", "scale+temp"]:
            for sf in ([1e-18] if case == "Conventional" else [1e-17, 1e-18, 1e-19]):
                f_, g_, b_ = mc_verify(A, sig_lev, sec, layouts[case], sf, mod, n_mc, rng, d_mm)
                rows.append({"case": case, "model": mod or "random", "sig_f": sf, "formal_rmse_mm": f_,
                             "mc_gauss_mm": g_, "mc_boot_mm": b_, "ratio_gauss": g_ / f_, "ratio_boot": b_ / f_})
    return pd.DataFrame(rows)