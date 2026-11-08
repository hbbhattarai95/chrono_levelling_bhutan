"""Monte Carlo verification of the exact error propagation."""
import numpy as np
import pandas as pd
from .adjustment import error_cov, summarise


def mc_verify(A, sig_lev, sec, nodes, sig_f, model, n_mc, rng, d_mm):
    """Return (formal RMSE, Gaussian MC RMSE, bootstrap MC RMSE) in mm.

    Bootstrap: the random levelling part is resampled from the observed normalised
    misclosures (d_k / 2 sqrt(L_k)) sqrt(L_s) with random sign; the rest is Gaussian."""
    Cx, G, S = error_cov(A, sig_lev, model, nodes, sig_f, "aware", sec)
    formal = summarise(Cx)[0]["rmse_mm"]
    eg = rng.standard_normal((n_mc, len(S))) @ np.linalg.cholesky(S).T
    gauss = 1e3 * np.sqrt(np.mean((eg @ G.T) ** 2))
    m = A.shape[0]
    S_rest = S.copy(); S_rest[:m, :m] -= np.diag(sig_lev ** 2)
    w, V = np.linalg.eigh(S_rest); w = np.clip(w, 0, None)
    er = rng.standard_normal((n_mc, len(S))) @ (V * np.sqrt(w)).T
    zk = d_mm * 1e-3 / (2 * np.sqrt(sec.L_km.values))
    er[:, :m] += rng.choice(zk, (n_mc, m)) * rng.choice([-1.0, 1.0], (n_mc, m)) * np.sqrt(sec.L_km.values)
    boot = 1e3 * np.sqrt(np.mean((er @ G.T) ** 2))
    return formal, gauss, boot


def run_verification(A, sig_lev, sec, layouts, d_mm, n_mc=5000, rng=None):
    rng = rng or np.random.default_rng(20261104)
    rows = []
    for case in ["Conventional", "Case 1", "Case 2", "Case 3"]:
        for mod in ["", "scale+temp"]:
            for sf in ([1e-18] if case == "Conventional" else [1e-17, 1e-18, 1e-19]):
                f_, g_, b_ = mc_verify(A, sig_lev, sec, layouts[case], sf, mod, n_mc, rng, d_mm)
                rows.append({"case": case, "model": mod or "random", "sig_f": sf, "formal_rmse_mm": f_,
                             "mc_gauss_mm": g_, "mc_boot_mm": b_, "ratio_gauss": g_ / f_, "ratio_boot": b_ / f_})
    return pd.DataFrame(rows)
