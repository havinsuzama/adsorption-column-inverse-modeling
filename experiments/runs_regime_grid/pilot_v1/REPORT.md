# 256 → 512-cell refits for two regimes

Fast-kinetics and high-affinity regimes, first concentration-noise draw and first porosity observation, both designs. Each fit was restarted from its 256-cell estimate and re-optimized on a 512-cell grid (at most 15 objective evaluations). The data were generated on 1024 cells, and the 150 mg/L predictions are evaluated on 1024 cells for both parameter sets. These cases were picked after the regime pilot, so this is a grid check rather than an independent replicate.

| Regime | Design | k shift | P shift | KL shift | eps shift | 150 mg/L RMS: 256 → 512 |
|---|---|---:|---:|---:|---:|---:|
| fast_kinetics | one_condition | +0.243% | +0.356% | -2.192% | +0.001% | 0.0191987 → 0.0181488 |
| fast_kinetics | two_concentrations | -0.797% | +0.026% | -0.088% | +0.014% | 0.0043836 → 0.0047976 |
| high_affinity | one_condition | +0.114% | +0.122% | -2.029% | -0.000% | 0.0047328 → 0.0049749 |
| high_affinity | two_concentrations | -0.117% | +0.019% | -0.319% | +0.004% | 0.0058042 → 0.0059064 |

The shifts are changes between the 256- and 512-cell estimates, not errors against the true parameters. True-parameter errors, objectives, stopping status and run times are in result.json.

## Notes

- The design ranking did not change: two concentrations stay better at fast kinetics (0.0181 vs 0.0048) and one concentration stays better at high affinity (0.0050 vs 0.0059).
- The largest shift is in KL for the one-condition design (about −2.2% and −2.0%); all two-condition shifts are below 1%.
- All four fits stopped successfully with no active bounds. This covers a single noise pairing, so it says nothing about convergence across all realizations.
