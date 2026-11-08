# Simulation-Based Evaluation of Levelling Errors and Optical-Clock Networks for the Height System of Bhutan

Code and results for the manuscript

> Bhattarai, H. B. *Optical clocks against random and systematic levelling errors: clock spacing and siting for an 822 km Himalayan levelling line in Bhutan.* Submitted to *Journal of Geodesy*.

Repository: <https://github.com/hbbhattarai95/chrono_levelling_bhutan>

## Overview

This project investigates whether a small number of optical clocks can improve the accuracy and reliability of a conventional levelling line in Bhutan. It uses the digital levelling route from **Samdrup Jongkhar to Phuntsholing** (316 benchmarks, 315 independent sections, 822.0 km, heights 169–3734 m) and evaluates optical-clock configurations by simulation and generalised least-squares (GLS) adjustment.

Clocks do not replace levelling. Clock-derived height differences are treated as additional observations, and the study asks how they change the uncertainty of the adjusted benchmark heights.

Workflow:

1. Prepare the levelling network and benchmark data.
2. Estimate the random error of the levelling from forward–backward misclosures.
3. Recover the rod-temperature correction applied in the field and build systematic error models.
4. Simulate optical-clock comparisons with prescribed fractional-frequency uncertainty.
5. Combine levelling and clock observations in a GLS adjustment and propagate errors exactly.
6. Compare clock configurations, error models and clock uncertainties.
7. Derive a clock-spacing rule and search for the best clock locations.
8. Verify the formal uncertainties by Gaussian and bootstrap Monte Carlo simulation.
9. Export tables and figures.

## 1. Input data

| File | Description |
| --- | --- |
| `data/Digital_Levelling_Samdrup_Jongkhar_to_Phuntsholing.xlsx` | Section length, temperature, forward/backward height differences, mean, temperature-calibrated mean, tolerance |
| `data/Coordinates_for_BMs.csv` | Benchmark names, Bhutan National Grid coordinates, published heights |

The levelling data were provided by the National Land Commission Secretariat (NLCS) of Bhutan. If the data files are not included in this repository, they are available from the author on reasonable request with the permission of NLCS.

The notebook rebuilds the network in benchmark order with **Samdrup Jongkhar BM** as the fixed datum. The return leg of the Changangkha spur is the exact negative of the outbound leg (a copy, not an independent observation) and is removed (2.56 km). Four control points listed after the line are not part of it and are excluded. An optional Bhutan boundary file (`.geojson` or `.shp`) can be set in the notebook for the map in Fig. 1.

## 2. Levelling errors

### Random error

For section $i$ the forward–backward misclosure is $d_i = F_i + B_i$. The random error of the mean of the two runs is

$$
\delta = \sqrt{\frac{1}{4n}\sum_{i=1}^{n}\frac{d_i^2}{L_i}} \quad [\mathrm{mm}/\sqrt{\mathrm{km}}],
$$

with $L_i$ the section length in km. A bootstrap over sections gives the 95 % interval. All sections satisfy the 12 mm$\sqrt{\mathrm{km}}$ acceptance tolerance, so a normal distribution truncated at that limit is also fitted (sensitivity case). Each section is assigned $\sigma_{\mathrm{lev},i} = \delta\sqrt{L_i}$.

Result: $\delta$ = 2.53 mm/√km (2.35–2.70); truncated-normal estimate 2.78 mm/√km.

### Rod-temperature correction

The *Calibrated Mean* in the workbook equals the mean corrected by $10.0\ \mathrm{ppm/^\circ C}\times(T-10.0\ ^\circ\mathrm{C})$. The notebook recovers this coefficient by regression and accumulates the applied correction along the line: +203 mm at km 543 and −255 mm at Phuntsholing BM. An error in this correction propagates as a systematic error.

## 3. Optical-clock observation model

$$
\frac{\Delta f}{f} \approx \frac{\Delta W}{c^2} \approx \frac{g_0\,\Delta H}{c^2},
\qquad
\sigma_h = \frac{c^2}{g_0}\,\sigma_f ,
$$

so $\sigma_f = 10^{-18}$ corresponds to $\sigma_h$ = 9.16 mm. Clock uncertainties tested: $\sigma_f$ = 1×10⁻¹⁷, 3×10⁻¹⁸, 1×10⁻¹⁸, 3×10⁻¹⁹, 1×10⁻¹⁹.

All comparisons are made against a reference clock at Samdrup Jongkhar BM, so they share its error:

$$
\boldsymbol{\Sigma}_c = \sigma_h^2\left(\mathbf{I} + \mathbf{1}\mathbf{1}^{\mathrm{T}}\right).
$$

## 4. Clock configurations

| Configuration | Clock locations |
| --- | --- |
| Conventional | none (levelling only) |
| Case 1 | Samdrup Jongkhar BM, Phuntsholing BM |
| Case 2 | Case 1 + Changangkha Zero Order Point (km 667, 2437 m) |
| Case 3 | Case 1 + Thrumshingla Pass (km 361, 3734 m) |

## 5. Combined adjustment

$$
\mathbf{l} = \mathbf{A}\mathbf{x} + \boldsymbol{\varepsilon},
\qquad
\hat{\mathbf{x}} = \mathbf{G}\,\mathbf{l},
\quad
\mathbf{G} = \left(\mathbf{A}^{\mathrm{T}}\mathbf{W}\mathbf{A}\right)^{-1}\mathbf{A}^{\mathrm{T}}\mathbf{W}.
$$

Because the estimator is linear, the error covariance of the adjusted heights under the true covariance $\boldsymbol{\Sigma}$ is exactly

$$
\mathbf{C}_{\hat{x}} = \mathbf{G}\,\boldsymbol{\Sigma}\,\mathbf{G}^{\mathrm{T}} .
$$

- **Aware** estimator: $\mathbf{W} = \boldsymbol{\Sigma}^{-1}$ (systematic errors modelled), so $\mathbf{C}_{\hat{x}} = (\mathbf{A}^{\mathrm{T}}\boldsymbol{\Sigma}^{-1}\mathbf{A})^{-1}$.
- **Naive** estimator: weights from the random model only, as in routine adjustments; its true error is still computed with the full $\boldsymbol{\Sigma}$.

Levelling and clock observations are mutually independent; the shared reference clock is in $\boldsymbol{\Sigma}_c$. Accuracy measures: RMSE over all benchmarks, standard deviation at Phuntsholing BM, worst benchmark standard deviation.

## 6. Error models

| Model | Content |
| --- | --- |
| `random` | $\delta\sqrt{L_i}$ only |
| `scale` | + rod-scale error common to the line, $\sigma_\kappa$ = 10 ppm |
| `temp` | + 10 % error of the applied temperature coefficient |
| `scale+temp` | reference systematic scenario (both) |
| `block` | + independent 10 ppm scale error per 100 km block |

## 7. Sensitivity analyses

Levelling error coefficient (2.53 vs 2.78 mm/√km), clock uncertainty (10⁻¹⁷ to 10⁻¹⁹), rod-scale uncertainty (0–30 ppm), error model, estimator (naive vs aware), and number and location of clocks.

## 8. Clock-spacing rule

Between two clocks a distance $D$ apart, the worst-case standard deviation is approximately

$$
\sigma_{\max}^2 \approx \frac{\delta^2 D}{4} + \sigma_h^2 ,
\qquad
D \le \frac{4\left(\sigma_t^2-\sigma_h^2\right)}{\delta^2}
$$

for a target $\sigma_t$. The rule is checked against the full GLS for 2–17 equally spaced clocks (agreement within 1.2 mm) and inverted for targets of 10, 15, 20 and 25 mm.

## 9. Clock siting

- **Third-clock scan:** with clocks at both terminal benchmarks, a third clock is placed at each intermediate benchmark in turn (random and `scale+temp` models, $\sigma_f = 10^{-18}$).
- **Greedy selection:** starting from the two terminal clocks, clocks are added one at a time up to six, each time choosing the site that minimises the RMSE (`scale+temp`, aware, $\sigma_f = 10^{-18}$), from (i) all benchmarks and (ii) 13 road-accessible candidate sites (district towns and passes, listed in `functions/constants.py`).

Selected sites are computational candidates; deployment also needs access, power, a stable environment and a frequency link.

## 10. Monte Carlo verification

5000 realisations per scenario for 20 scenarios: Gaussian errors drawn from $\boldsymbol{\Sigma}$, and a bootstrap in which the random levelling part is resampled from the observed misclosures. Ratios of Monte Carlo to formal RMSE: 0.990–1.019 (Gaussian), 0.983–1.017 (bootstrap).

## 11. Key results ($\sigma_f = 10^{-18}$)

| | Conventional | Case 1 | Case 2 | Case 3 |
| --- | --- | --- | --- | --- |
| Random errors: RMSE / worst (mm) | 51.5 / 72.3 | 30.9 / 36.7 | 26.4 / 33.2 | 23.3 / 29.2 |
| + scale and temperature: RMSE / worst (mm) | 56.6 / 76.7 | 41.8 / 58.9 | 28.6 / 39.7 | 24.9 / 32.7 |

Clocks at Samdrup Jongkhar BM and Phuntsholing BM (similar heights) cannot detect height-proportional errors; one extra clock at a high pass makes the result almost insensitive to rod-scale error.

## 12. Outputs (`results/`)

| File | Contents |
| --- | --- |
| `scenario_results.csv` | All cases × error models × clock uncertainties × estimators |
| `verification.csv` | Monte Carlo verification |
| `greedy_siting.csv` | Greedy clock selection |
| `sensitivity_scale.csv` | Sensitivity to rod-scale uncertainty |
| `sensitivity_delta.csv` | Sensitivity to the levelling error coefficient |
| `spacing_rule_check.csv` | Spacing rule vs full GLS |
| `third_clock_scan.csv` | RMSE for every third-clock location |
| `calibration_summary.csv` | Levelling-error and temperature-model calibration |
| `stations.csv` | Benchmarks, coordinates, heights, chainage, cumulative temperature correction |
| `key_numbers.json` | Main numbers quoted in the manuscript |
| `figures/Fig1–Fig6` | 600 dpi PNG and vector PDF |

Figures: (1) study area, height profile and cumulative temperature correction; (2) height standard deviation along the line; (3) RMSE vs clock uncertainty and the spacing rule; (4) systematic errors and rod-scale sensitivity; (5) third-clock scan and greedy siting; (6) misclosure distribution and Monte Carlo verification.

## 13. Project structure

```text
chrono_levelling/
├── chrono_levelling.ipynb     main notebook (runs everything, about 1 minute)
├── data/                      input data
├── functions/
│   ├── constants.py           constants, clock uncertainties, systematic settings, candidate sites
│   ├── data.py                reading data, spur removal, chainage, temperature correction
│   ├── calibration.py         delta (bootstrap, truncated normal), temperature model
│   ├── adjustment.py          design matrix, covariances, GLS, exact error propagation
│   ├── scenarios.py           scenario runs and sensitivity analyses
│   ├── siting.py              third-clock scan, greedy selection, spacing rule
│   ├── verification.py        Gaussian and bootstrap Monte Carlo
│   ├── mapping.py             Bhutan National Grid projection, boundary reading
│   └── figures.py             Fig. 1–6
├── results/                   written by the notebook
├── requirements.txt
└── README.md
```

## 14. Running the analysis

```bash
git clone https://github.com/hbbhattarai95/chrono_levelling_bhutan.git
cd chrono_levelling_bhutan
pip install -r requirements.txt     # numpy, scipy, pandas, openpyxl, matplotlib, jupyter
jupyter notebook chrono_levelling.ipynb
```

Run all cells in order. `N_MC = 5000` and `SEED = 20261104` (first code cell) reproduce all numbers in the manuscript.

## 15. Limitations

- Clock observations are simulated; clock performance is represented by a prescribed fractional-frequency uncertainty that includes link, tidal and local-tie errors.
- Published heights are taken as truth, so results describe error propagation, not the accuracy of the existing heights.
- Gravity is not part of the error model.
- Systematic error magnitudes are assumptions, tested over 0–30 ppm.
- Levelling and clock observations are treated as independent; clock errors as white.
- Practical deployment needs engineering and site-feasibility studies.

## 16. Research question

**Can a small number of optical-clock observations improve the precision of a conventional levelling line in Bhutan, and how do clock uncertainty, levelling error characteristics and clock locations influence that improvement?**

## Citation

If you use this code, please cite the manuscript above (details will be updated on publication).
