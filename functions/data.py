"""Reading the levelling workbook and benchmark coordinates."""
import numpy as np
import pandas as pd


def load_data(xlsx, csv):
    """Return (stations, extra control points, sections, removed spur length in km).

    sections columns: i, j (station indices), L_km, T_C (temperature), F, B (forward /
    backward height differences), M (mean), C (temperature-calibrated mean), tol (tolerance),
    tcorr = C - M (applied temperature correction, m).
    """
    sec = pd.read_excel(xlsx, header=None, skiprows=5).iloc[:316]
    sta = pd.read_csv(csv)
    extra = sta[sta.Slno > 316].reset_index(drop=True)            # control points off the line
    sta = sta[sta.Slno <= 316].reset_index(drop=True)
    assert len(sec) == 316 and len(sta) == 316
    # The line visits Changlingmithang -> Changangkha -> Changlingmithang (out-and-back spur)
    path = np.array(list(range(0, 257)) + [258, 257, 258] + list(range(259, 316)))
    H = sta.Levelling_Height.values.astype(float)
    assert np.abs(sec[11].astype(float).values - H[path[1:]]).max() < 0.002
    f = lambda c: sec[c].astype(float).values
    s = pd.DataFrame({"i": path[:-1], "j": path[1:], "L_km": f(1), "T_C": f(0), "F": f(4), "B": f(5),
                      "M": f(6), "C": f(7), "tol": f(13)})
    # Return leg of the spur is the exact negative of the outbound leg: a copy -> removed
    assert abs(f(6)[257] + f(6)[258]) < 1e-9
    spur_km = float(s.loc[258, "L_km"])
    s = s.drop(index=258).reset_index(drop=True)
    s["tcorr"] = s.C - s.M
    chain = np.full(316, np.nan); chain[0] = 0.0
    for r in s.itertuples():
        if np.isnan(chain[r.j]):
            chain[r.j] = chain[r.i] + r.L_km
    sta["chain_km"] = chain
    return sta, extra, s, spur_km


def cumulative_temperature_correction(sta, sec):
    """Cumulative applied temperature correction at each station (mm)."""
    cum = np.zeros(len(sta)); seen = {0}
    for r in sec.itertuples():
        if r.j not in seen:
            cum[r.j] = cum[r.i] + r.tcorr; seen.add(r.j)
    return cum * 1e3


def find_station(sta, name):
    return int(np.where(sta.Name_BM.str.strip() == name)[0][0])
