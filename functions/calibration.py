"""Noise calibration from forward/backward misclosures and the temperature model."""
import numpy as np
from scipy import stats, optimize


def estimate_delta(sec, n_boot=5000, seed=1):
    """Random levelling error delta (mm/sqrt(km)) of the mean of forward and backward runs.

    delta^2 = mean(d^2 / L) / 4 with d = F + B.
    Returns delta, 95% bootstrap interval, misclosures d (mm), truncated-normal delta,
    tolerance (mm/sqrt(km)).
    """
    d = (sec.F.values + sec.B.values) * 1e3
    q = d ** 2 / sec.L_km.values
    delta = np.sqrt(q.mean() / 4)
    rng = np.random.default_rng(seed)
    bs = np.array([np.sqrt(rng.choice(q, len(q)).mean() / 4) for _ in range(n_boot)])
    # all sections lie within the tolerance -> fit a normal truncated at +-tol
    z = d / np.sqrt(sec.L_km.values)
    t = float((sec.tol.values / np.sqrt(sec.L_km.values)).mean() * 1e3)
    nll = lambda s_: -np.sum(stats.norm.logpdf(z, 0, s_) - np.log(2 * stats.norm.cdf(t, 0, s_) - 1))
    s_tr = optimize.minimize_scalar(nll, bounds=(1, 30), method="bounded").x
    return delta, np.percentile(bs, [2.5, 97.5]), d, s_tr / 2, t


def temperature_model_check(sec):
    """Fit ppm = a*T + b to the applied calibration. Returns (a ppm/degC, T0 degC, max residual ppm)."""
    ok = np.abs(sec.M) > 5
    ppm = (sec.C / sec.M - 1)[ok] * 1e6
    a, b = np.polyfit(sec.T_C[ok], ppm, 1)
    resid = np.abs(ppm - (a * sec.T_C[ok] + b)).max()
    return a, -b / a, resid
