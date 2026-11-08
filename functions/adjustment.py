"""Combined levelling and clock adjustment using GLS and exact error propagation."""

import numpy as np
from .constants import C2_G0, SIG_KAPPA, SIG_R, BLOCK_KM


def design(sec, n_sta):
    """Construct the levelling design matrix with station 0 fixed as the datum.

    Each observation represents a height difference:

        dH_ij = H_j - H_i

    Therefore, each row contains -1 for the starting station and +1
    for the ending station. The datum column is removed because its
    height is fixed to zero.
    """
    A = np.zeros((len(sec), n_sta))
    A[np.arange(len(sec)), sec.j.values] += 1
    A[np.arange(len(sec)), sec.i.values] -= 1

    return A[:, 1:]


def sys_columns(sec, model, sig_kappa=SIG_KAPPA,
                sig_r=SIG_R, block_km=BLOCK_KM):
    """Construct the design matrix for systematic levelling errors.

    Systematic errors are represented using a low-rank covariance model:

        Sigma_sys = U P U^T

    where U contains the sensitivity of each observation to the
    systematic error parameters, and P contains their variances.

    Supported error models:

    1. scale:
       A common scale-factor error affecting all levelling observations.

    2. temp:
       A temperature-related error represented by the temperature
       correction term.

    3. block:
       Independent scale-factor errors assigned to consecutive
       distance blocks of specified length.

    The total covariance is obtained by adding the systematic
    covariance to the random observation covariance.
    """
    dH = sec.C.values
    cols, var = [], []

    # Common rod-scale error affecting the entire levelling network.
    if "scale" in model:
        cols.append(dH)
        var.append(sig_kappa ** 2)

    # Temperature-related systematic error.
    if "temp" in model:
        cols.append(sec.tcorr.values)
        var.append(sig_r ** 2)

    # Independent rod-scale errors for successive distance blocks.
    if "block" in model:
        start = np.r_[0, np.cumsum(sec.L_km.values)[:-1]]
        b = (start // block_km).astype(int)

        for k in np.unique(b):
            cols.append(np.where(b == k, dH, 0.0))
            var.append(sig_kappa ** 2)

    # Return empty arrays when only random errors are considered.
    if not cols:
        return np.zeros((len(sec), 0)), np.zeros(0)

    return np.column_stack(cols), np.array(var)


def clock_block(nodes, nunk, sig_h):
    """Construct the design matrix and covariance for clock comparisons.

    All clock comparisons are made relative to the reference clock
    located at nodes[0].

    Each observation represents a height difference between the
    reference clock and another clock.

    The covariance accounts for the uncertainty of both clocks.
    Since the same reference clock is used in every comparison,
    its uncertainty introduces correlations between observations:

        Sigma_c = sigma_h^2 (I + 11^T)

    where:
        sigma_h : standard uncertainty of an individual clock
        I       : identity matrix
        11^T    : common reference-clock contribution
    """
    Kc = len(nodes) - 1
    Ac = np.zeros((Kc, nunk))

    # Each row represents a comparison with the reference clock.
    for r, nd in enumerate(nodes[1:]):
        Ac[r, nd - 1] = 1.0

    # Include the shared uncertainty of the reference clock.
    Sc = sig_h ** 2 * (np.eye(Kc) + np.ones((Kc, Kc)))

    return Ac, Sc


def full_system(A, sig_lev, U, p, nodes, sig_h):
    """Combine levelling and clock observations into one adjustment system.

    The levelling covariance consists of random and systematic
    contributions:

        Sigma_l = Sigma_random + U P U^T

    The combined observation vector contains both levelling
    observations and clock comparisons.

    Two covariance matrices are returned:

    Sigma:
        Complete covariance, including systematic levelling errors.

    Sigma_n:
        Random-only levelling covariance, retained for the naive
        estimator that ignores systematic error correlations.
    """
    # Full levelling observation covariance.
    Sl = np.diag(sig_lev ** 2) + (U * p) @ U.T

    # If no clock comparisons are available, use levelling alone.
    if nodes is None or len(nodes) < 2:
        return A, Sl, np.diag(sig_lev ** 2)

    # Construct the clock observation system.
    Ac, Sc = clock_block(nodes, A.shape[1], sig_h)

    m, k = A.shape[0], Ac.shape[0]

    # Stack levelling and clock design matrices.
    Ab = np.vstack([A, Ac])

    # Combined covariance, assuming levelling and clock errors
    # are mutually independent.
    S = np.zeros((m + k, m + k))
    S[:m, :m] = Sl
    S[m:, m:] = Sc

    # Alternative covariance that ignores systematic levelling errors.
    Sn = np.zeros_like(S)
    Sn[:m, :m] = np.diag(sig_lev ** 2)
    Sn[m:, m:] = Sc

    return Ab, S, Sn


def gls(Ab, Sw):
    """Calculate the GLS gain matrix and parameter covariance.

    The generalised least-squares solution is:

        x_hat = (A^T Sigma^-1 A)^-1 A^T Sigma^-1 l

    The gain matrix is:

        G = (A^T W A)^-1 A^T W

    and the covariance of the estimated parameters, assuming Sw
    is the true observation covariance, is:

        Q_x = (A^T W A)^-1

    where W = Sw^-1.
    """
    W = np.linalg.inv(Sw)

    # Inverse of the normal equation matrix.
    Ninv = np.linalg.inv(Ab.T @ W @ Ab)

    # Matrix mapping observation errors into parameter errors.
    G = Ninv @ Ab.T @ W

    return G, Ninv


def error_cov(A, sig_lev, model, nodes, sig_f,
              estimator="aware", sec=None, **kw):
    """Propagate the complete observation covariance into adjusted heights.

    The adjusted parameter errors are related to observation errors by:

        delta_x = G delta_l

    Therefore, the exact covariance of the adjusted heights is:

        C_x = G Sigma G^T

    Two estimation strategies can be investigated:

    aware:
        Uses the complete covariance, including systematic
        levelling errors, when calculating the GLS weights.

    naive:
        Uses only random levelling errors when calculating the
        weights, but propagates the complete covariance when
        evaluating the resulting height uncertainties.

    The clock uncertainty is converted from fractional frequency
    uncertainty into an equivalent height uncertainty using:

        sigma_h = c^2 / g_0 * sigma_f

    Returns:
        Cx     : covariance matrix of adjusted heights
        G      : GLS gain matrix
        Sigma  : complete observation covariance matrix
    """
    # Construct systematic error sensitivity matrix and variances.
    U, p = sys_columns(sec, model, **kw)

    # Build the combined levelling and clock system.
    Ab, S, Sn = full_system(
        A, sig_lev, U, p, nodes, C2_G0 * sig_f
    )

    # Select the covariance used to calculate the estimator weights.
    Sw = S if estimator == "aware" else Sn

    # Calculate the gain matrix using the selected estimator.
    G, _ = gls(Ab, Sw)

    # Propagate the complete observation covariance.
    Cx = G @ S @ G.T

    return Cx, G, S


def summarise(Cx):
    """Summarise adjusted-height uncertainties in millimetres.

    The diagonal elements of Cx represent the variances of
    individual adjusted station heights.

    Standard deviations are calculated as:

        sigma_i = sqrt(Cx_ii)

    The reported RMSE is the root-mean-square of station
    standard deviations:

        RMSE = sqrt(mean(sigma_i^2))

    The datum station is included in the returned standard
    deviation profile with zero uncertainty.
    """
    # Convert standard deviations from metres to millimetres.
    sd = 1e3 * np.sqrt(np.clip(np.diag(Cx), 0, None))

    # Calculate network-wide uncertainty statistics.
    summary = {
        "rmse_mm": float(np.sqrt(np.mean(sd ** 2))),
        "end_mm": float(sd[-1]),
        "worst_mm": float(sd.max())
    }

    # Include the fixed datum station at the beginning.
    profile = np.r_[0, sd]

    return summary, profile