# Sparse KL profile at high affinity

High-affinity regime, first concentration-noise draw and first porosity observation, 512-cell fits. KL was fixed at the fitted value plus 0, ±1 and ±2 local standard deviations (clipped to the original bounds), and k, P and eps were re-optimized at each point. Five points per design, with a second start at the two ends: 14 fits. Points and starts were written to protocol.json before fitting.

With δ = log(KL) − log(KL_hat) and s the local log(KL) standard deviation, a quadratic objective would rise by (δ/s)². The table compares this with the observed rise 2ΔJ relative to the re-optimized centre.

| Design | Requested z | Fixed KL | Actual z | 2ΔJ observed | Quadratic prediction | Starts | Start gap |
|---|---:|---:|---:|---:|---:|---:|---:|
| one_condition | -2 | 0.055704 | -2.0000 | 4.88276 | 4.00000 | 2 | 2.85e-10 |
| one_condition | -1 | 0.093550 | -1.0000 | 1.07868 | 1.00000 | 1 | 0 |
| one_condition | +0 | 0.157110 | +0.0000 | 0.00000 | 0.00000 | 1 | 0 |
| one_condition | +1 | 0.263854 | +1.0000 | 0.75482 | 1.00000 | 1 | 0 |
| one_condition | +2 | 0.400000 | +1.8025 | 2.06384 | 3.24909 | 2 | 2.49e-08 |
| two_concentrations | -2 | 0.117482 | -2.0000 | 5.30612 | 4.00000 | 2 | 1.04e-10 |
| two_concentrations | -1 | 0.138102 | -1.0000 | 1.16824 | 1.00000 | 1 | 0 |
| two_concentrations | +0 | 0.162342 | +0.0000 | 0.00000 | 0.00000 | 1 | 0 |
| two_concentrations | +1 | 0.190837 | +1.0000 | 0.89238 | 1.00000 | 1 | 0 |
| two_concentrations | +2 | 0.224333 | +2.0000 | 3.10438 | 4.00000 | 2 | 6.19e-10 |

## Notes

- In both designs the objective rises faster than the quadratic prediction toward low KL and more slowly toward high KL, so the local standard deviation is not a symmetric interval here.
- When a point hits the KL upper bound (0.4), its actual z is below the requested value and the quadratic prediction uses the actual z.
- Points with an objective below the centre (beyond 1e-6): 0. Unsuccessful stops: 0/14; numerical exceptions: 0.
- Five points per design show the shape of the profile; they are not a confidence interval.
