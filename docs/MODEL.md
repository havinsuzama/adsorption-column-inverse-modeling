# Model, units and defaults

```
C_t + v*C_z - D*C_zz + (1000*rho_b/eps)*q_t = 0
q_t = k*(qmax*KL*C/(1+KL*C) - q)
```

`C`: mg/L pore liquid; `q`: mg/g dry adsorbent; `rho_b`: kg dry adsorbent/L bed; `eps`: porosity; `D`: m²/min; `k`: min⁻¹; `KL`: L/mg; `qmax`: mg/g. The factor 1000 converts kg to g. `v = Q/(area·eps)`, with Q converted from mL/min to m³/min.

Initial condition C = q = 0. Inlet: `C - (D/v)*C_z = C0` (Danckwerts); outlet: `C_z = 0`.

Default geometry: length 0.06 m, diameter 0.01 m, Q = 0.157 mL/min, duration 24000 min. Synthetic defaults: eps = 0.4, rho_b = 0.8 kg/L, D = 7.496197819628419e-6 m²/min, k = 0.001 min⁻¹, qmax = 60 mg/g, KL = 0.04 L/mg. qref = 48.1 mg/g is only a normalization scale. Time is normalized by the run duration.

## Numerics

Finite volumes with upwind advection, centered dispersion and BDF time integration with a sparse Jacobian. The outlet concentration is the last cell average. Grid accuracy is checked by comparing 1024- and 2048-cell solutions.

## Identifiability

Scaling rho_b → a·rho_b, qmax → qmax/a, q → q/a leaves all liquid concentrations unchanged. The inverse problems therefore estimate log(k, P, KL, eps) with P = rho_b·qmax, and keep density fixed. Dispersion is kept fixed when the flow changes.
