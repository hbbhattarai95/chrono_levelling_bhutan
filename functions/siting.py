"""Clock siting (third-clock scan, greedy selection) and the clock-spacing rule."""
import numpy as np
import pandas as pd
from .constants import C2_G0
from .adjustment import error_cov, summarise


def scan_third_clock(A, sig_lev, sec, n, sig_f, model):
    """RMSE with clocks at both terminal benchmarks and a third clock at each benchmark k."""
    out = np.full(n, np.nan)
    for k in range(1, n - 1):
        Cx, _, _ = error_cov(A, sig_lev, model, np.array([0, k, n - 1]), sig_f, "aware", sec)
        out[k] = summarise(Cx)[0]["rmse_mm"]
    return out


def greedy(A, sig_lev, sec, candidates, base, sig_f, model, kmax, metric="rmse_mm"):
    """Add clocks one at a time, each time choosing the candidate that minimises `metric`."""
    chosen, hist = list(base), []
    Cx, _, _ = error_cov(A, sig_lev, model, np.array(chosen), sig_f, "aware", sec)
    hist.append((list(chosen), summarise(Cx)[0]))
    for _ in range(kmax):
        best = None
        for c in candidates:
            if c in chosen:
                continue
            nodes = np.array(sorted(chosen[1:] + [c])); nodes = np.r_[0, nodes[nodes != 0]]
            Cx, _, _ = error_cov(A, sig_lev, model, nodes, sig_f, "aware", sec)
            s = summarise(Cx)[0]
            if best is None or s[metric] < best[1][metric]:
                best = (c, s)
        chosen.append(best[0]); hist.append((list(chosen), best[1]))
    return hist


def greedy_table(sta, hist, label):
    rows = []
    for nodes, s in hist:
        names = [f"{sta.Name_BM[x].strip()} (km {sta.chain_km[x]:.0f})" for x in sorted(nodes, key=lambda z: sta.chain_km[z])]
        rows.append({"candidates": label, "n_clocks": len(nodes), "sites": "; ".join(names), **s})
    return pd.DataFrame(rows)


def spacing_rule(delta_mm, D_km, sig_h_mm):
    """Worst-case s.d. between two clocks D km apart: sqrt(delta^2 D / 4 + sig_h^2)."""
    return np.sqrt(delta_mm ** 2 * D_km / 4 + sig_h_mm ** 2)


def spacing_check(A, sig_lev, sec, sta, delta_mm, sig_f, nclks=(2, 3, 4, 5, 7, 9, 13, 17)):
    """Compare the rule with the full adjustment for equally spaced clocks."""
    ch = sta.chain_km.values; rows = []
    for k in nclks:
        nodes = np.unique([int(np.argmin(np.abs(ch - t))) for t in np.linspace(0, ch.max(), k)])
        Cx, _, _ = error_cov(A, sig_lev, "", nodes, sig_f, "aware", sec)
        D = float(np.diff(ch[nodes]).max())
        rows.append({"sig_f": sig_f, "n_clocks": len(nodes), "D_km": D, "gls_worst_mm": summarise(Cx)[0]["worst_mm"],
                     "rule_mm": float(spacing_rule(delta_mm, D, C2_G0 * sig_f * 1e3))})
    return pd.DataFrame(rows)
