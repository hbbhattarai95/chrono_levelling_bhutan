"""Physical constants and default settings for the study."""

# Fundamental physical constants.
C_LIGHT = 299_792_458.0                  # Speed of light in vacuum (m/s)
G0 = 9.80665                             # Conventional gravitational acceleration (m/s²)

# Conversion factor between fractional frequency difference and height difference.
# From the relativistic relation: ΔW / c² ≈ Δf / f
# Using ΔW ≈ g₀ Δh gives: Δh = (c² / g₀) × (Δf / f)
C2_G0 = C_LIGHT ** 2 / G0                # Height conversion factor (m)

# Fractional frequency uncertainties considered for the optical clocks.
# These values represent different clock performance levels.
SIGMA_F = [
    1e-17,
    3e-18,
    1e-18,
    3e-19,
    1e-19
]

# Uncertainty assigned to the levelling instrument's rod-scale calibration.
SIG_KAPPA = 10e-6                        # 1-sigma scale uncertainty (10 ppm)

# Relative uncertainty of the applied temperature correction.
SIG_R = 0.10                             # 10% uncertainty

# Length of each section used in the block-wise scale-error model.
BLOCK_KM = 100.0                         # Block length (km)

# Candidate locations for optical clock installation.
# Only benchmarks with suitable road access are included.
# Names follow the benchmark identifiers in Coordinates_for_BMs.csv.
TOWNS = [
    "Samdrup Jongkhar BM",
    "Kheri BM",
    "Mongar BM",
    "Lingmithang BM",
    "Thrumshingla Pass",
    "Ura Junction",
    "B?M near Chamkhar chhu bridge",
    "B?M near Trongsa Dzong",
    "Pelela Pass BM",
    "Wangdue BM",
    "Dochula Pass BM",
    "Changangkha Zero Order Point",
    "Chukha Dam BM",
    "Gedu BM",
    "Phuntsholing BM"
]