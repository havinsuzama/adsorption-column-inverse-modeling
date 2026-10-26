# Four parameter regimes

Each regime changes one parameter from the base case (k = 0.001, P = 48, KL = 0.04, eps = 0.4): slow kinetics k = 0.0005, fast kinetics k = 0.002, low affinity KL = 0.01, high affinity KL = 0.16. Per regime: 2 concentration-noise draws × 2 porosity observations × 2 designs × 2 starts = 16 fits; 64 in total. The same standard noise draws are used across regimes and designs, so the four cells per regime are not independent replicates. Both designs have 40 observations.

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

- slow_kinetics: two concentrations better in 2/4 pairs for P, 2/4 for KL and 2/4 for the 150 mg/L prediction. Unsuccessful stops 0/16, active bounds 0/16, largest start gap 1.73e-09.
- fast_kinetics: two concentrations better in 4/4 pairs for P, 4/4 for KL and 4/4 for the 150 mg/L prediction. Unsuccessful stops 0/16, active bounds 0/16, largest start gap 8.43e-09.
- low_affinity: two concentrations better in 2/4 pairs for P, 2/4 for KL and 4/4 for the 150 mg/L prediction. Unsuccessful stops 0/16, active bounds 0/16, largest start gap 3.18e-09.
- high_affinity: two concentrations better in 2/4 pairs for P, 3/4 for KL and 2/4 for the 150 mg/L prediction. Unsuccessful stops 0/16, active bounds 0/16, largest start gap 2.91e-08.

## Notes

- In each cell the start with the lower observed objective was kept, never the one closer to the true parameters. All starts, including unsuccessful ones, are in each regime's result.json; paired differences are in comparison.json.
- Reference data use 1024 cells and were checked against 2048 cells for 100, 150 and 200 mg/L before any fit started (reference_checks.json in each regime folder). Fits use 256 cells.
- The 150 mg/L prediction is computed on the 1024-cell grid. This condition was used during development, so it is not a blind test.
- P and porosity are the same in every regime, so their effect is not tested. qmax is not identified separately from bulk density.
