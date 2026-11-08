"""Physical constants and default study settings."""
C_LIGHT = 299_792_458.0                  # m/s
G0 = 9.80665                             # m/s^2, conventional gravity
C2_G0 = C_LIGHT ** 2 / G0                # m per unit fractional frequency (9.165e15 m)

SIGMA_F = [1e-17, 3e-18, 1e-18, 3e-19, 1e-19]   # clock fractional uncertainties studied
SIG_KAPPA = 10e-6                        # rod-scale uncertainty (1 sigma), 10 ppm
SIG_R = 0.10                             # relative error of the applied temperature correction
BLOCK_KM = 100.0                         # block length of the block-wise scale model

# Candidate clock sites with road access (benchmark names as in Coordinates_for_BMs.csv)
TOWNS = ["Samdrup Jongkhar BM", "Kheri BM", "Mongar BM", "Lingmithang BM", "Thrumshingla Pass",
         "Ura Junction", "B?M near Chamkhar chhu bridge", "B?M near Trongsa Dzong", "Pelela Pass BM",
         "Wangdue BM", "Dochula Pass BM", "Changangkha Zero Order Point", "Chukha Dam BM",
         "Gedu BM", "Phuntsholing BM"]
