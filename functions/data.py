"""Read levelling observations and benchmark coordinates."""

import numpy as np
import pandas as pd


def load_data(xlsx, csv):
    """Load the levelling workbook and benchmark coordinate data.

    The function reads the levelling observations, separates additional
    control points from benchmarks along the main levelling line, and
    reconstructs the observation sequence.

    The Changlingmithang–Changangkha spur is traversed in both directions.
    Since the return observation duplicates the outbound height difference,
    it is removed before the adjustment.

    The resulting section table contains:

        i, j     : indices of the starting and ending benchmarks
        L_km     : levelling section length (km)
        T_C      : observed temperature (degree C)
        F, B     : forward and backward height differences (m)
        M        : mean height difference (m)
        C        : temperature-calibrated mean (m)
        tol      : allowable levelling tolerance (m)
        tcorr    : applied temperature correction, C - M (m)

    Returns:
        sta      : benchmark information, including cumulative chainage
        extra    : additional control points outside the main line
        s        : processed levelling observations
        spur_km  : length of the removed return section (km)
    """

    # Read the levelling observations, skipping the workbook header.
    # Only the first 316 observations belong to the study line.
    sec = pd.read_excel(xlsx, header=None, skiprows=5).iloc[:316]

    # Load benchmark coordinates and associated station information.
    sta = pd.read_csv(csv)

    # Separate additional control points that are not part of the
    # main levelling sequence.
    extra = sta[sta.Slno > 316].reset_index(drop=True)
    sta = sta[sta.Slno <= 316].reset_index(drop=True)

    # Check that the number of observations and main-line benchmarks
    # agrees with the expected network configuration.
    assert len(sec) == 316 and len(sta) == 316

    # Reconstruct the actual levelling route.
    # The sequence passes through Changlingmithang, visits Changangkha,
    # and returns to Changlingmithang before continuing along the line.
    path = np.array(
        list(range(0, 257))
        + [258, 257, 258]
        + list(range(259, 316))
    )

    # Extract the benchmark levelling heights.
    H = sta.Levelling_Height.values.astype(float)

    # Verify that the workbook's observed heights match the benchmark
    # heights in the reconstructed station sequence.
    assert np.abs(
        sec[11].astype(float).values - H[path[1:]]
    ).max() < 0.002

    # Helper function to retrieve a workbook column as numerical data.
    f = lambda c: sec[c].astype(float).values

    # Assemble the levelling observations into a structured table.
    s = pd.DataFrame({
        "i": path[:-1],
        "j": path[1:],
        "L_km": f(1),
        "T_C": f(0),
        "F": f(4),
        "B": f(5),
        "M": f(6),
        "C": f(7),
        "tol": f(13)
    })

    # Confirm that the return observation of the spur is the exact
    # negative of the corresponding outbound observation.
    assert abs(f(6)[257] + f(6)[258]) < 1e-9

    # Record the length of the duplicated return section before removing it.
    spur_km = float(s.loc[258, "L_km"])

    # Remove the redundant return observation and reset the row indices.
    s = s.drop(index=258).reset_index(drop=True)

    # Calculate the temperature correction applied to each section.
    s["tcorr"] = s.C - s.M

    # Calculate cumulative chainage from the first benchmark.
    # The starting benchmark is assigned zero chainage, and each
    # subsequent station inherits the chainage of its predecessor
    # plus the length of the connecting levelling section.
    chain = np.full(316, np.nan)
    chain[0] = 0.0

    for r in s.itertuples():
        if np.isnan(chain[r.j]):
            chain[r.j] = chain[r.i] + r.L_km

    # Store the cumulative distance along the levelling line.
    sta["chain_km"] = chain

    return sta, extra, s, spur_km


def cumulative_temperature_correction(sta, sec):
    """Calculate the cumulative temperature correction at each benchmark.

    The correction for a station is obtained by accumulating the
    applied section corrections along the levelling line:

        TC_j = TC_i + tcorr_ij

    The first benchmark is used as the reference and assigned
    zero cumulative correction.

    Returns:
        Cumulative temperature correction at each station (mm).
    """

    # Initialise all station corrections to zero.
    cum = np.zeros(len(sta))

    # Track stations whose cumulative correction has been calculated.
    seen = {0}

    # Accumulate section corrections along the levelling sequence.
    for r in sec.itertuples():
        if r.j not in seen:
            cum[r.j] = cum[r.i] + r.tcorr
            seen.add(r.j)

    # Convert the cumulative corrections from metres to millimetres.
    return cum * 1e3


def find_station(sta, name):
    """Return the index of a benchmark identified by its name.

    Leading and trailing whitespace in benchmark names is ignored.
    """
    return int(
        np.where(sta.Name_BM.str.strip() == name)[0][0]
    )