"""Journal-style figures (Times-compatible serif, 600 dpi PNG + vector PDF).

Every figure function takes the results dictionary R built in the notebook and an output
folder, saves PNG + PDF and returns the matplotlib Figure (so it also shows in the notebook).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from .constants import C2_G0, SIGMA_F
from .siting import spacing_rule

COL = {
    "conv": "#000000",
    "c1": "#0072B2",
    "c2": "#D55E00",
    "c3": "#009E73",
    "clock": "#CC79A7",
    "grey": "#888888",
}
LS = {"conv": "-", "c1": "--", "c2": "-.", "c3": ":"}
W2 = 7.2  # double-column width (inches)
CASE_LABEL = {
    "c1": "Case 1: Samdrup Jongkhar BM + Phuntsholing BM",
    "c2": "Case 2: Case 1 + Changangkha",
    "c3": "Case 3: Case 1 + Thrumshingla",
}


def set_style():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Liberation Serif", "Times New Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 8,
            "axes.labelsize": 8.5,
            "axes.titlesize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.2,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linewidth": 0.4,
            "lines.linewidth": 1.1,
            "savefig.dpi": 600,
            "pdf.fonttype": 42,
            "figure.dpi": 110,
        }
    )


def save(fig, base):
    fig.savefig(str(base) + ".png", bbox_inches="tight", pad_inches=0.03)
    fig.savefig(str(base) + ".pdf", bbox_inches="tight", pad_inches=0.03)


def panel(ax, s):
    ax.text(
        -0.02,
        1.03,
        s,
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=9.5,
        va="bottom",
        ha="right",
    )


def _sites(R):
    n = len(R["sta"])
    return 0, R["br"], R["thr"], n - 1


def fig1_study_area(R, F):
    sta, extra = R["sta"], R["extra"]
    ch = sta.chain_km.values
    hgt = sta.Levelling_Height.values
    n = len(sta)
    s0, br, th, sn = _sites(R)
    names = {s0: "Samdrup Jongkhar BM", br: "Changangkha", sn: "Phuntsholing BM"}
    fig = plt.figure(figsize=(W2, 5.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.2, 1], hspace=0.30)
    ax = fig.add_subplot(gs[0])
    ax.grid(False)
    if R.get("boundary"):
        for r in R["boundary"]:
            ax.fill(
                r[:, 0] / 1e3,
                r[:, 1] / 1e3,
                color="#eeeeee",
                ec="#777777",
                lw=0.7,
                zorder=0,
            )
    ax.plot(sta.Easting / 1e3, sta.Northing / 1e3, color="#9a9a9a", lw=0.6, zorder=1)
    ax.scatter(
        sta.Easting / 1e3,
        sta.Northing / 1e3,
        s=4,
        color="k",
        zorder=2,
        label=f"Benchmarks (n = {n})",
    )
    ax.scatter(
        extra.Easting / 1e3,
        extra.Northing / 1e3,
        s=10,
        marker="s",
        color=COL["c3"],
        zorder=2,
        label=f"Control points off the line (n = {len(extra)})",
    )
    tw = [t for t in R["towns_idx"] if t not in (s0, br, sn)]
    ax.scatter(
        sta.Easting.values[tw] / 1e3,
        sta.Northing.values[tw] / 1e3,
        s=22,
        marker="o",
        facecolor="none",
        ec=COL["c1"],
        lw=0.8,
        zorder=3,
        label="Candidate clock sites",
    )
    ax.scatter(
        sta.Easting.values[[s0, br, sn]] / 1e3,
        sta.Northing.values[[s0, br, sn]] / 1e3,
        s=60,
        marker="*",
        color=COL["c2"],
        ec="k",
        lw=0.4,
        zorder=4,
        label="Clock sites, Cases 1-2",
    )
    for s_, (dx, dy, ha) in {
        s0: (-8, -2, "right"),
        br: (-8, 0, "right"),
        sn: (6, -10, "left"),
    }.items():
        ax.annotate(
            names[s_],
            (sta.Easting.values[s_] / 1e3, sta.Northing.values[s_] / 1e3),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7.2,
            ha=ha,
        )
    ax.annotate(
        "Thrumshingla",
        (sta.Easting.values[th] / 1e3, sta.Northing.values[th] / 1e3),
        xytext=(0, 12),
        textcoords="offset points",
        fontsize=7.2,
        ha="center",
        bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.8),
    )
    ax.set_aspect("equal")
    ax.set_xlabel("Easting (km)")
    ax.set_ylabel("Northing (km)")
    ax.set_xlim(150, 420)
    ax.legend(frameon=False, loc="center", bbox_to_anchor=(0.6, 0.42), fontsize=7)
    panel(ax, "(a)")
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(ch, hgt, color="k", lw=0.9, label="Height")
    ax2.scatter(
        ch[[s0, br, sn]],
        hgt[[s0, br, sn]],
        s=45,
        marker="*",
        color=COL["c2"],
        ec="k",
        lw=0.4,
        zorder=3,
    )
    for s_, (dx, dy, ha) in {
        s0: (8, 4, "left"),
        br: (8, 6, "left"),
        sn: (-6, 6, "right"),
    }.items():
        ax2.annotate(
            names[s_],
            (ch[s_], hgt[s_]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7.2,
            ha=ha,
        )
    ax2.annotate(
        "Thrumshingla 3734 m",
        (ch[th], hgt[th]),
        xytext=(6, 2),
        textcoords="offset points",
        fontsize=7.2,
    )
    ax2.set_xlabel("Distance along the line (km)")
    ax2.set_ylabel("Height (m)")
    panel(ax2, "(b)")
    ax3 = ax2.twinx()
    ax3.spines["right"].set_visible(True)
    ax3.grid(False)
    ax3.plot(
        ch,
        R["cum_corr_mm"],
        color=COL["c1"],
        lw=0.9,
        ls="--",
        label="Cumulative temperature correction",
    )
    ax3.set_ylabel("Cumulative temperature\ncorrection (mm)", color=COL["c1"])
    ax3.tick_params(axis="y", colors=COL["c1"])
    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax3.get_legend_handles_labels()
    ax2.legend(h1 + h2, l1 + l2, frameon=False, loc="upper right", fontsize=7)
    save(fig, F / "Fig1_study_area")
    return fig


def fig2_error_along_line(R, F):
    ch = R["sta"].chain_km.values
    prof = R["prof"]
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.8), sharey=True)
    for a, sf, lab in zip(axs, [1e-18, 1e-19], ["(a)", "(b)"]):
        a.plot(
            ch,
            prof[("conv", "", 1e-17)],
            color=COL["conv"],
            ls=LS["conv"],
            label="Conventional levelling",
        )
        for c in ["c1", "c2", "c3"]:
            a.plot(
                ch,
                prof[(c, "", sf)],
                color=COL[c],
                ls=LS[c],
                lw=1.3 if c == "c3" else 1.1,
                label=CASE_LABEL[c],
            )
        for s_ in _sites(R):
            a.axvline(ch[s_], color=COL["clock"], lw=0.5, alpha=0.6)
        a.set_xlabel("Distance along the line (km)")
        a.set_title(r"$\sigma_f = 10^{%d}$" % int(np.log10(sf)))
        panel(a, lab)
    axs[0].set_ylabel("Height standard deviation (mm)")
    axs[0].set_ylim(0, 80)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(
        h,
        l,
        frameon=False,
        loc="lower center",
        ncol=2,
        bbox_to_anchor=(0.5, -0.17),
        fontsize=7.2,
    )
    save(fig, F / "Fig2_error_along_line")
    return fig


def fig3_clock_uncertainty_and_spacing(R, F):
    res = R["res"]
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.8))
    a = axs[0]
    a.axhline(R["base_rmse"], color="k", label="Conventional levelling")
    for c, k in [("Case 1", "c1"), ("Case 2", "c2"), ("Case 3", "c3")]:
        d = res[(res.case == c) & (res.model == "random")].sort_values("sig_f")
        a.plot(d.sig_f, d.rmse_mm, color=COL[k], ls=LS[k], marker="o", ms=3, label=c)
    a.set_xscale("log")
    a.invert_xaxis()
    a.set_ylim(0, R["base_rmse"] * 1.12)
    a.set_xlabel(r"Clock fractional uncertainty $\sigma_f$")
    a.set_ylabel("Height RMSE, all benchmarks (mm)")
    a.legend(frameon=False, loc="lower left")
    panel(a, "(a)")
    sec_ax = a.secondary_xaxis(
        "top", functions=(lambda x: x * C2_G0 * 1e3, lambda y: y / (C2_G0 * 1e3))
    )
    sec_ax.set_xlabel(r"$\sigma_h$ of one clock (mm)", fontsize=7.5)
    a = axs[1]
    D = np.linspace(5, 900, 300)
    sp = R["spacing"]
    for sf, k in [(1e-18, "c2"), (1e-19, "c1")]:
        a.plot(
            D,
            spacing_rule(R["delta"], D, C2_G0 * sf * 1e3),
            color=COL[k],
            ls=LS[k],
            label=r"Eq. (7), $\sigma_f = 10^{%d}$" % int(np.log10(sf)),
        )
        d = sp[sp.sig_f == sf]
        a.scatter(d.D_km, d.gls_worst_mm, s=14, color=COL[k], ec="k", lw=0.3, zorder=3)
    for y in (10, 20):
        a.axhline(y, color=COL["grey"], lw=0.6, ls=":")
    a.set_xlabel("Clock spacing $D$ (km)")
    a.set_ylabel("Worst-case height s.d. (mm)")
    a.set_xlim(0, 900)
    a.set_ylim(0, 40)
    a.legend(frameon=False, loc="lower right")
    panel(a, "(b)")
    fig.tight_layout()
    save(fig, F / "Fig3_clock_uncertainty_and_spacing")
    return fig


def fig4_systematic_errors(R, F):
    ch = R["sta"].chain_km.values
    prof, ps, sv = R["prof"], R["prof_sys"], R["sens"]
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.9))
    a = axs[0]
    a.plot(
        ch,
        prof[("conv", "scale+temp", 1e-17)],
        color="k",
        label="Conventional levelling",
    )
    a.plot(
        ch,
        prof[("conv", "", 1e-17)],
        color="k",
        lw=0.6,
        ls=":",
        label="Conventional, random errors",
    )
    for c in ["c1", "c2", "c3"]:
        a.plot(
            ch,
            ps[(c, "aware")],
            color=COL[c],
            ls=LS[c],
            lw=1.3 if c == "c3" else 1.1,
            label=CASE_LABEL[c],
        )
    for s_ in _sites(R):
        a.axvline(ch[s_], color=COL["clock"], lw=0.5, alpha=0.6)
    a.set_xlabel("Distance along the line (km)")
    a.set_ylabel("Height standard deviation (mm)")
    a.set_ylim(0, 85)
    a.legend(
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.2),
        ncol=2,
        fontsize=6.6,
    )
    panel(a, "(a)")
    a = axs[1]
    for case, k in [
        ("Conventional", "conv"),
        ("Case 1", "c1"),
        ("Case 2", "c2"),
        ("Case 3", "c3"),
    ]:
        d = sv[sv.case == case]
        a.plot(
            d.sig_kappa_ppm,
            d.worst_mm,
            color=COL[k],
            ls=LS[k],
            marker="o",
            ms=3,
            label=case,
        )
    a.set_xlabel(r"Rod-scale uncertainty $\sigma_\kappa$ (ppm)")
    a.set_ylabel("Worst benchmark s.d. (mm)")
    a.set_ylim(0, None)
    a.legend(frameon=False, loc="lower right", ncol=2)
    panel(a, "(b)")
    fig.tight_layout()
    save(fig, F / "Fig4_systematic_errors")
    return fig


def fig5_clock_siting(R, F):
    sta = R["sta"]
    ch = sta.chain_km.values
    hgt = sta.Levelling_Height.values
    _, br, th, _ = _sites(R)
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.8))
    a = axs[0]
    a.plot(ch, R["scan"]["random"], color=COL["c1"], label="Random errors")
    a.plot(
        ch,
        R["scan"]["scale+temp"],
        color=COL["c2"],
        ls="--",
        label="Random + systematic",
    )
    for key, k in [("random", "c1"), ("scale+temp", "c2")]:
        i = int(np.nanargmin(R["scan"][key]))
        a.scatter(
            ch[i], R["scan"][key][i], color=COL[k], s=20, zorder=3, ec="k", lw=0.3
        )
    a.axvline(ch[br], color=COL["clock"], lw=0.6)
    a.axvline(ch[th], color=COL["clock"], lw=0.6, ls="--")
    top = np.nanmax(R["scan"]["scale+temp"])
    a.text(ch[br] + 6, top, "Changangkha", fontsize=6.8, color=COL["clock"])
    a.text(
        ch[th] - 6, top, "Thrumshingla", fontsize=6.8, color=COL["clock"], ha="right"
    )
    a.set_xlabel("Chainage of the third clock (km)")
    a.set_ylabel("Height RMSE (mm)")
    ax2 = a.twinx()
    ax2.fill_between(ch, hgt, color="#dddddd", zorder=0, lw=0)
    ax2.set_ylim(0, 12000)
    ax2.axis("off")
    a.set_zorder(ax2.get_zorder() + 1)
    a.patch.set_visible(False)
    a.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 0.1))
    panel(a, "(a)")
    a = axs[1]
    for key, k, lab in [
        ("all", "c1", "Any benchmark"),
        ("towns", "c2", "Candidate sites"),
    ]:
        hst = R["greedy"][key]
        x = range(2, 2 + len(hst))
        a.plot(
            x,
            [h[1]["rmse_mm"] for h in hst],
            color=COL[k],
            marker="o",
            ms=3,
            label=lab + ", RMSE",
        )
        a.plot(
            x,
            [h[1]["worst_mm"] for h in hst],
            color=COL[k],
            marker="s",
            ms=3,
            ls="--",
            label=lab + ", worst",
        )
    a.set_xlabel("Number of clocks")
    a.set_ylabel("Height error (mm)")
    a.set_ylim(0, None)
    a.set_xticks(range(2, 7))
    a.legend(frameon=False, loc="upper right", fontsize=6.6)
    panel(a, "(b)")
    fig.tight_layout()
    save(fig, F / "Fig5_clock_siting")
    return fig


def fig6_model_checks(R, F):
    sec = R["sec"]
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.7))
    a = axs[0]
    z = R["d_mm"] / np.sqrt(sec.L_km.values)
    a.hist(z, bins=30, color="#bbbbbb", ec="k", lw=0.4, density=True)
    xg = np.linspace(-16, 16, 300)
    t = R["tol_mm"]
    for s_, k, lab, trunc in [
        (2 * R["delta"], "c2", r"Normal, $\delta$ = %.2f" % R["delta"], False),
        (
            2 * R["delta_tr"],
            "c1",
            r"Truncated normal, $\delta$ = %.2f" % R["delta_tr"],
            True,
        ),
    ]:
        pdf = stats.norm.pdf(xg, 0, s_)
        if trunc:
            pdf = np.where(np.abs(xg) <= t, pdf / (2 * stats.norm.cdf(t, 0, s_) - 1), 0)
        a.plot(xg, pdf, color=COL[k], lw=1.1, label=lab + r" mm/$\sqrt{\mathrm{km}}$")
    for v in (-t, t):
        a.axvline(v, color=COL["grey"], lw=0.6, ls=":")
    a.set_xlabel(r"$(F+B)/\sqrt{L}$ (mm/$\sqrt{\mathrm{km}}$)")
    a.set_ylabel("Density")
    a.set_xlim(-16, 16)
    a.set_ylim(0, 0.14)
    a.legend(frameon=False, loc="upper left", fontsize=6.6)
    panel(a, "(a)")
    a = axs[1]
    v = R["verif"]
    lim = max(v.formal_rmse_mm.max(), v.mc_gauss_mm.max(), v.mc_boot_mm.max()) * 1.05
    a.plot([0, lim], [0, lim], color="k", lw=0.7)
    a.scatter(
        v.formal_rmse_mm,
        v.mc_gauss_mm,
        s=18,
        color=COL["c1"],
        ec="k",
        lw=0.3,
        label="Gaussian Monte Carlo",
        zorder=3,
    )
    a.scatter(
        v.formal_rmse_mm,
        v.mc_boot_mm,
        s=18,
        marker="^",
        color=COL["c2"],
        ec="k",
        lw=0.3,
        label="Bootstrap Monte Carlo",
        zorder=3,
    )
    a.set_xlabel("Formal RMSE (mm)")
    a.set_ylabel("Monte Carlo RMSE (mm)")
    a.set_aspect("equal")
    a.set_xlim(0, lim)
    a.set_ylim(0, lim)
    a.legend(frameon=False, loc="upper left")
    panel(a, "(b)")
    fig.tight_layout()
    save(fig, F / "Fig6_model_checks")
    return fig
