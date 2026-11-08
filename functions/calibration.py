"""Calibration of levelling noise and validation of the temperature model."""

import numpy as np
from scipy import stats, optimize


def estimate_delta(sec, n_boot=5000, seed=1):
    """Estimate random levelling error from forward and backward observations.

    The random error of the mean of the forward and backward runs is
    estimated from their misclosures:

        d = F + B

    The corresponding random-error coefficient is:

        delta^2 = (1 / 4) * mean(d^2 / L)

    where:
        d      : forward-backward misclosure (m)
        L      : levelling section length (km)
        delta  : random error coefficient (mm / sqrt(km))

    A bootstrap procedure is used to estimate the 95% confidence
    interval of delta by repeatedly resampling the observed sections.

    Since all observed sections satisfy the prescribed tolerance,
    a second estimate is obtained by fitting a zero-mean normal
    distribution truncated at the section-specific tolerance limits.

    Returns:
        delta       : estimated random error coefficient
        ci          : 95% bootstrap confidence interval
        d           : forward-backward misclosures (mm)
        delta_tr    : truncated-normal estimate
        t           : mean normalised tolerance (mm / sqrt(km))
    """

    # Calculate forward-backward misclosures and convert them to mm.
    d = (sec.F.values + sec.B.values) * 1e3

    # Normalise the squared misclosures by their section lengths.
    q = d ** 2 / sec.L_km.values

    # Estimate the random error coefficient from the mean misclosure.
    # The factor 1/4 accounts for averaging the two independent runs.
    delta = np.sqrt(q.mean() / 4)

    # Bootstrap the observed sections to quantify estimation uncertainty.
    rng = np.random.default_rng(seed)

    bs = np.array([
        np.sqrt(rng.choice(q, len(q)).mean() / 4)
        for _ in range(n_boot)
    ])

    # Normalise each misclosure by the square root of its section length.
    z = d / np.sqrt(sec.L_km.values)

    # Calculate the mean normalised tolerance in mm / sqrt(km).
    t = float(
        (sec.tol.values / np.sqrt(sec.L_km.values)).mean() * 1e3
    )

    # Fit a zero-mean normal distribution truncated at +/- tolerance.
    # The likelihood is conditional on the observations falling
    # within the prescribed tolerance limits.
    nll = lambda s_: -np.sum(
        stats.norm.logpdf(z, 0, s_)
        - np.log(2 * stats.norm.cdf(t, 0, s_) - 1)
    )

    # Estimate the standard deviation by minimising the negative
    # log-likelihood over the specified search interval.
    s_tr = optimize.minimize_scalar(
        nll,
        bounds=(1, 30),
        method="bounded"
    ).x

    # Return the estimates, bootstrap interval and misclosures.
    return (
        delta,
        np.percentile(bs, [2.5, 97.5]),
        d,
        s_tr / 2,
        t
    )


def temperature_model_check(sec):
    """Check the applied rod-temperature calibration using levelling data.

    The applied calibration is examined by fitting a linear relationship:

        ppm = a * T + b

    where:
        ppm : relative calibration correction in parts per million
        T   : observed temperature (degree C)
        a   : temperature sensitivity (ppm / degree C)
        b   : intercept (ppm)

    The reference temperature at which the fitted correction becomes
    zero is:

        T0 = -b / a

    The maximum absolute residual indicates how closely the applied
    calibration follows the fitted linear temperature relationship.

    Returns:
        a       : fitted temperature coefficient (ppm / degree C)
        T0      : estimated reference temperature (degree C)
        resid   : maximum absolute fitting residual (ppm)
    """

    # Exclude observations with very small applied corrections,
    # as these can make the relative calibration ratio unstable.
    ok = np.abs(sec.M) > 5

    # Express the applied correction relative to the measured value
    # as a parts-per-million quantity.
    ppm = (sec.C / sec.M - 1)[ok] * 1e6

    # Fit a straight line between temperature and calibration correction.
    a, b = np.polyfit(sec.T_C[ok], ppm, 1)

    # Calculate the largest absolute departure from the fitted model.
    resid = np.abs(
        ppm - (a * sec.T_C[ok] + b)
    ).max()

    # Estimate the reference temperature from the zero-crossing.
    return a, -b / a, resid