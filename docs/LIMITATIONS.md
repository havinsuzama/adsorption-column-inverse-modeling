# Limitations

- All data are synthetic. There are no measured radionuclide breakthrough curves and no isotope-specific validation.
- The model is single-species Langmuir–LDF. Competing ions, pH effects, radioactive decay and detector response are not included.
- The four regimes vary only kinetics and affinity; capacity and porosity are the same in all of them.
- Crossed cells share noise draws, so they are not independent replicates. Medians are descriptive and no confidence intervals are reported.
- Every fit used two fixed starts. All fits stopped on the `ftol`/`xtol` criteria; two starts agreeing does not prove a global optimum.
- The 150 mg/L condition was looked at during development, so it is not a blind test.
- Reference data and fits come from the same finite-volume model at different resolutions. This tests the numerics, not the physics.
- Local (Gauss–Newton) standard deviations ignore nonlinearity and parameter bounds. The KL profile shows that the nonlinearity matters.
- The porosity-precision calculation keeps the fitted state fixed. It is not a simulated new measurement followed by a refit.
- The PINN has not passed its forward-solution check.
