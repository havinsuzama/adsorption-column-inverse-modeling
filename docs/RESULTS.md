# Results

All numbers come from the archived synthetic runs in `experiments/`. Medians are descriptive. Crossed cells share noise draws, so they are not independent replicates. Errors are relative to the true parameters used to generate the data.

## 1. Crossed noise pilot (36 fits)

One true parameter set (k = 0.001 min⁻¹, P = 48, KL = 0.04 L/mg, ε = 0.4). Three concentration-noise draws were crossed with three auxiliary porosity observations. Two designs with the same number of points: 40 points at 100 mg/L, or 20 + 20 points at 100 and 200 mg/L. Two optimizer starts per fit; the start with the lower objective was kept.

| Design | Median abs. k error | P error | KL error | ε error | Median 150 mg/L outlet RMS |
|---|---:|---:|---:|---:|---:|
| one_condition | 10.95% | 5.27% | 35.45% | 4.69% | 0.016172 |
| two_concentrations | 4.36% | 1.89% | 9.18% | 4.45% | 0.003481 |

Two concentrations gave lower absolute errors for all four parameters and a lower 150 mg/L prediction error in 9 of 9 pairs. Details: [crossed pilot report](../experiments/runs_crossed_replication/pilot_v1/REPORT.md).

## 2. Four parameter regimes (64 fits)

Each regime changes one parameter from the base case: slow kinetics (k = 0.0005), fast kinetics (k = 0.002), low affinity (KL = 0.01), high affinity (KL = 0.16). Per regime: 2 noise draws × 2 porosity observations × 2 designs × 2 starts.

| Regime | Design | Median abs. k error | P error | KL error | Median 150 mg/L outlet RMS |
|---|---|---:|---:|---:|---:|
| slow_kinetics | one_condition | 10.16% | 4.33% | 38.84% | 0.012344 |
| slow_kinetics | two_concentrations | 1.80% | 1.86% | 6.47% | 0.005758 |
| fast_kinetics | one_condition | 20.71% | 7.53% | 67.68% | 0.028928 |
| fast_kinetics | two_concentrations | 6.55% | 0.52% | 2.37% | 0.004083 |
| low_affinity | one_condition | 11.75% | 10.54% | 28.77% | 0.012528 |
| low_affinity | two_concentrations | 2.82% | 2.86% | 6.10% | 0.005263 |
| high_affinity | one_condition | 7.54% | 0.69% | 14.24% | 0.006001 |
| high_affinity | two_concentrations | 6.80% | 0.32% | 3.87% | 0.005889 |

The two-concentration prediction was better in 2/4 pairs for slow kinetics, 4/4 for fast kinetics, 4/4 for low affinity and 2/4 for high affinity. At high affinity the median improvement is only about 1.9%. So the benefit of a second concentration depends on the regime.

Reference solutions use 1024 cells and were checked against 2048 cells; fits use 256 cells. Details: [regime report](../experiments/runs_parameter_regimes/pilot_v1/REPORT.md).

## 3. Grid check: 256 → 512 cells

For the first noise pairing of the fast-kinetics and high-affinity regimes, each fit was restarted from its 256-cell estimate on a 512-cell grid.

| Regime | Design | 150 mg/L outlet RMS after refit |
|---|---|---:|
| fast_kinetics | one_condition | 0.0181488 |
| fast_kinetics | two_concentrations | 0.0047976 |
| high_affinity | one_condition | 0.0049749 |
| high_affinity | two_concentrations | 0.0059064 |

The ranking of the designs did not change. The largest shift was in KL for the one-condition design (about 2%). Details: [grid report](../experiments/runs_regime_grid/pilot_v1/REPORT.md).

## 4. Local uncertainty, KL profile and porosity precision

- **Local uncertainty.** At the four 512-cell fits, two concentrations reduced the local log(KL) standard deviation from 0.433 to 0.070 (fast kinetics) and from 0.518 to 0.162 (high affinity). The log(ε) standard deviation stayed at about 0.05, the assumed precision of the porosity observation. [Report](../experiments/runs_local_regime_uncertainty/pilot_v1/REPORT.md)
- **KL profile.** At high affinity, KL was fixed at ±1 and ±2 local standard deviations and the other parameters were re-optimized (14 fits). The objective rises faster than the quadratic prediction toward low KL (4.88 vs 4.00 at z = −2, one condition) and slower toward high KL (2.06 vs 3.25 at the bound-clipped point). [Report](../experiments/runs_high_affinity_profile/pilot_v1/REPORT.md)
- **Porosity precision.** Keeping the fitted states fixed, tightening the porosity measurement from σ = 0.05 to 0.01 (log scale) reduced the local log(KL) standard deviation by at most 0.16%. [Report](../experiments/runs_auxiliary_precision_geometry/pilot_v1/REPORT.md)

## 5. PINN

Only an early-time PINN diagnostic is included. It has not passed its forward-solution check, and nothing here compares it with the finite-volume solver.

## Reproducing

`run_demo.py` uses smaller grids and a smaller optimization budget than the experiments, so it will not reproduce these numbers. To rerun the experiments, see [experiments/README.md](../experiments/README.md).
