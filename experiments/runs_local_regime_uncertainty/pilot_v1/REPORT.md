# Local uncertainty at four selected 512-cell fits

For each of the four 512-cell fits (fast kinetics and high affinity, both designs) I computed the Jacobian of the outlet curve with respect to x = log(k, P, KL, eps), divided by the noise standard deviation, plus one row [0, 0, 0, 20] for the porosity observation (log-sigma 0.05). The local covariance is (JᵀJ)⁻¹, computed with an SVD. Noise scales are the known synthetic values and were not rescaled from the residuals. This is the Gauss–Newton approximation at the fitted point, not a nonlinear confidence region.

| Regime | Design | log(k) std | log(P) std | log(KL) std | log(eps) std |
|---|---|---:|---:|---:|---:|
| fast_kinetics | one_condition | 0.21224 | 0.06826 | 0.43254 | 0.05000 |
| fast_kinetics | two_concentrations | 0.09673 | 0.01196 | 0.07044 | 0.04989 |
| high_affinity | one_condition | 0.11033 | 0.02977 | 0.51845 | 0.05000 |
| high_affinity | two_concentrations | 0.04601 | 0.00904 | 0.16171 | 0.04999 |

## Notes

- Two concentrations reduce the local log(KL) standard deviation about six-fold at fast kinetics (0.433 → 0.070) and about three-fold at high affinity (0.518 → 0.162).
- log(eps) std stays at about 0.05, the precision of the porosity observation, so almost all porosity information comes from that observation.
- Without the porosity row (same fitted state, no refit), the one-condition problem is poorly conditioned: condition numbers 5250 (fast kinetics) and 7652 (high affinity), against 211 and 251 with the row. Singular values, correlations and these concentration-only results are in result.json.
- In the high-affinity pair, the two-concentration fit still has the larger 150 mg/L outlet RMS (0.0059 vs 0.0050). Local uncertainty and the error of one noise realization measure different things.
- The two designs are evaluated at different fitted points, so the table does not isolate the pure design effect. Small log-std values can be read roughly as relative uncertainty; large ones cannot.
