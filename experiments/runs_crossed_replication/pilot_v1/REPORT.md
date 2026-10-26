# Crossed noise pilot

One true parameter set (k = 0.001 min⁻¹, P = 48, KL = 0.04 L/mg, eps = 0.4). Three concentration-noise draws were crossed with three auxiliary porosity observations. Two designs with 40 points each (40 at 100 mg/L, or 20 + 20 at 100 and 200 mg/L) and two optimizer starts: 36 fits. In each cell the start with the lower observed objective was kept. The nine cells share noise draws, so they are not nine independent datasets.

| Design | Median abs. k error | P error | KL error | eps error | Median 150 mg/L outlet RMS |
|---|---:|---:|---:|---:|---:|
| one_condition | 10.95% | 5.27% | 35.45% | 4.69% | 0.016172 |
| two_concentrations | 4.36% | 1.89% | 9.18% | 4.45% | 0.003481 |

Paired comparison, absolute error with two concentrations minus one:

- k_min: lower with two concentrations in 9/9 pairs; median paired difference -6.60 percentage points.
- density_capacity_product: lower with two concentrations in 9/9 pairs; median paired difference -3.38 percentage points.
- KL_l_mg: lower with two concentrations in 9/9 pairs; median paired difference -22.75 percentage points.
- porosity: lower with two concentrations in 9/9 pairs; median paired difference -0.24 percentage points.
- 150 mg/L outlet RMS: lower with two concentrations in 9/9 pairs.

## Notes

- The porosity observations were [0.4186, 0.4284, 0.418]. All three happen to lie above 0.4; no bias was added to the measurement model.
- Reference data use 1024 cells and were checked against 2048 cells; fits use 256 cells. All three feeds passed the preset checks (outlet RMS < 0.001, max difference < 0.002, relative mass error < 1e-6); see reference_checks.json.
- The 150 mg/L prediction is computed with the fitted parameters on the 1024-cell grid. This condition was used during development, so it is not a blind test.
- 0 of 36 fits reported an unsuccessful stop and 0 hit a parameter bound. The largest objective gap between the two starts was 1.09e-08.
- Density was fixed at 0.8254 kg/L, so P = rho_b·qmax is the capacity result; qmax alone is not identified.
- Mean errors per factor level and interaction residuals are in factor_effects.json. Three levels per factor are too few for variance components.
