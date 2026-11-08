"""Combined levelling + clock adjustment (generalised least squares) and exact error propagation."""
import numpy as np
from .constants import C2_G0, SIG_KAPPA, SIG_R, BLOCK_KM


def design(sec, n_sta):
    """Levelling design matrix; station 0 (Samdrup Jongkhar BM) is the datum."""
    A = np.zeros((len(sec), n_sta))
    A[np.arange(len(sec)), sec.j.values] += 1
    A[np.arange(len(sec)), sec.i.values] -= 1
    return A[:, 1:]


def sys_columns(sec, model, sig_kappa=SIG_KAPPA, sig_r=SIG_R, block_km=BLOCK_KM):
    """Systematic part of the levelling covariance: Sigma_sys = U diag(p) U^T.

    model is a string containing any of 'scale', 'temp', 'block' ('' = random only).
    """
    dH = sec.C.values; cols, var = [], []
    if "scale" in model:
        cols.append(dH); var.append(sig_kappa ** 2)
    if "temp" in model:
        cols.append(sec.tcorr.values); var.append(sig_r ** 2)
    if "block" in model:
        start = np.r_[0, np.cumsum(sec.L_km.values)[:-1]]
        b = (start // block_km).astype(int)
        for k in np.unique(b):
            cols.append(np.where(b == k, dH, 0.0)); var.append(sig_kappa ** 2)
    if not cols:
        return np.zeros((len(sec), 0)), np.zeros(0)
    return np.column_stack(cols), np.array(var)


def clock_block(nodes, nunk, sig_h):
    """Design and covariance of clock comparisons against the reference clock at nodes[0]."""
    Kc = len(nodes) - 1
    Ac = np.zeros((Kc, nunk))
    for r, nd in enumerate(nodes[1:]):
        Ac[r, nd - 1] = 1.0
    Sc = sig_h ** 2 * (np.eye(Kc) + np.ones((Kc, Kc)))     # common reference-clock error
    return Ac, Sc


def full_system(A, sig_lev, U, p, nodes, sig_h):
    """Stacked design, true covariance and random-only (naive) covariance."""
    Sl = np.diag(sig_lev ** 2) + (U * p) @ U.T
    if nodes is None or len(nodes) < 2:
        return A, Sl, np.diag(sig_lev ** 2)
    Ac, Sc = clock_block(nodes, A.shape[1], sig_h)
    m, k = A.shape[0], Ac.shape[0]
    Ab = np.vstack([A, Ac])
    S = np.zeros((m + k, m + k)); S[:m, :m] = Sl; S[m:, m:] = Sc
    Sn = np.zeros_like(S); Sn[:m, :m] = np.diag(sig_lev ** 2); Sn[m:, m:] = Sc
    return Ab, S, Sn


def gls(Ab, Sw):
    """GLS gain matrix G and (A^T W A)^-1 for weight covariance Sw."""
    W = np.linalg.inv(Sw)
    Ninv = np.linalg.inv(Ab.T @ W @ Ab)
    return Ninv @ Ab.T @ W, Ninv


def error_cov(A, sig_lev, model, nodes, sig_f, estimator="aware", sec=None, **kw):
    """Exact error covariance G Sigma G^T of the adjusted heights.

    estimator: 'aware' (weights from full covariance) or 'naive' (random model weights).
    Returns (Cx, G, Sigma).
    """
    U, p = sys_columns(sec, model, **kw)
    Ab, S, Sn = full_system(A, sig_lev, U, p, nodes, C2_G0 * sig_f)
    G, _ = gls(Ab, S if estimator == "aware" else Sn)
    return G @ S @ G.T, G, S


def summarise(Cx):
    """RMSE, end-point and worst standard deviation (mm), plus the s.d. profile incl. datum."""
    sd = 1e3 * np.sqrt(np.clip(np.diag(Cx), 0, None))
    return {"rmse_mm": float(np.sqrt(np.mean(sd ** 2))), "end_mm": float(sd[-1]),
            "worst_mm": float(sd.max())}, np.r_[0, sd]
