# Porosity-measurement precision at fixed fits

At the four 512-cell fits, each fitted state and its concentration Jacobian were kept fixed and only the assumed precision of the porosity observation was changed: I_new = I_base + (1/σ_new² − 1/0.05²) e_eps e_epsᵀ. The covariance was updated with a rank-one formula and checked against a direct solve (largest relative difference 5e-14). No new model solves, observations or fits were made.

| Regime | Design | Assumed log-eps σ | log(k) std | log(P) std | log(KL) std | log(eps) std |
|---|---|---:|---:|---:|---:|---:|
| fast_kinetics | one_condition | 0.01 | 0.21223 | 0.06814 | 0.43185 | 0.01000 |
| fast_kinetics | one_condition | 0.02 | 0.21223 | 0.06816 | 0.43194 | 0.02000 |
| fast_kinetics | one_condition | 0.05 | 0.21224 | 0.06826 | 0.43254 | 0.05000 |
| fast_kinetics | one_condition | 0.10 | 0.21227 | 0.06863 | 0.43469 | 0.10000 |
| fast_kinetics | one_condition | 0.20 | 0.21239 | 0.07007 | 0.44320 | 0.19997 |
| fast_kinetics | two_concentrations | 0.01 | 0.09636 | 0.01195 | 0.07043 | 0.01000 |
| fast_kinetics | two_concentrations | 0.02 | 0.09640 | 0.01195 | 0.07043 | 0.01999 |
| fast_kinetics | two_concentrations | 0.05 | 0.09673 | 0.01196 | 0.07044 | 0.04989 |
| fast_kinetics | two_concentrations | 0.10 | 0.09789 | 0.01199 | 0.07050 | 0.09914 |
| fast_kinetics | two_concentrations | 0.20 | 0.10211 | 0.01208 | 0.07069 | 0.19335 |
| high_affinity | one_condition | 0.01 | 0.11033 | 0.02974 | 0.51800 | 0.01000 |
| high_affinity | one_condition | 0.02 | 0.11033 | 0.02974 | 0.51806 | 0.02000 |
| high_affinity | one_condition | 0.05 | 0.11033 | 0.02977 | 0.51845 | 0.05000 |
| high_affinity | one_condition | 0.10 | 0.11034 | 0.02986 | 0.51986 | 0.10000 |
| high_affinity | one_condition | 0.20 | 0.11039 | 0.03024 | 0.52545 | 0.19998 |
| high_affinity | two_concentrations | 0.01 | 0.04600 | 0.00903 | 0.16163 | 0.01000 |
| high_affinity | two_concentrations | 0.02 | 0.04600 | 0.00903 | 0.16164 | 0.02000 |
| high_affinity | two_concentrations | 0.05 | 0.04601 | 0.00904 | 0.16171 | 0.04999 |
| high_affinity | two_concentrations | 0.10 | 0.04606 | 0.00906 | 0.16197 | 0.09993 |
| high_affinity | two_concentrations | 0.20 | 0.04625 | 0.00914 | 0.16299 | 0.19942 |

## Result

Going from σ = 0.05 to 0.01 reduces the local log(KL) standard deviation by:

- fast kinetics, one condition: 0.160%
- fast kinetics, two concentrations: 0.024%
- high affinity, one condition: 0.087%
- high affinity, two concentrations: 0.051%

A more precise porosity measurement mostly narrows porosity itself and does almost nothing for affinity at these points.

## Notes

- This is a local calculation at fixed fits. It does not simulate new measurements and refits, and it does not cover systematic measurement bias.
- The two designs sit at different fitted points, so differences between designs are not only due to information content.
- The KL profile is not quadratic, so these numbers should not be turned into confidence intervals. Density was fixed; qmax uncertainty was not computed.
