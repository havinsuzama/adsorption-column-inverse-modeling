# Experiment scripts and archived runs

Each `runs_*/pilot_v1/` folder holds one experiment: the settings fixed before the run (`protocol.json`), the synthetic observations, every fit including both optimizer starts, figures and a `REPORT.md`. Run the scripts from this folder.

| Order | Script | Output folder | What it does |
|---|---|---|---|
| 1 | `crossed_replication.py` | `runs_crossed_replication` | 3 × 3 crossed noise pilot, 36 fits |
| 2 | `parameter_regimes.py` | `runs_parameter_regimes` | Four parameter regimes, 64 fits |
| 3 | `refine_selected_regimes.py` | `runs_regime_grid` | 256 → 512-cell refits of two regimes |
| 4 | `local_regime_uncertainty.py` | `runs_local_regime_uncertainty` | Local (Gauss–Newton) uncertainty |
| 5 | `high_affinity_profile.py` | `runs_high_affinity_profile` | Sparse KL profile, 14 fits |
| 6 | `auxiliary_precision_geometry.py` | `runs_auxiliary_precision_geometry` | Effect of porosity-measurement precision |

`nonlinear_design_pilot.py` is an older, independent sensitivity implementation kept as a numerical cross-check.

## Regenerating the reports (no fitting)

```bash
python report_crossed_replication.py
python report_parameter_regimes.py
python report_high_affinity_profile.py
```

## Rerunning the experiments

The fitting scripts (1, 2, 3 and 5) can take hours. Each script stores SHA-256 hashes of its sources in `provenance.json` and refuses to resume if they have changed; `fits.json` and `partial.json` are the resume checkpoints. Copy a run folder before rerunning it, and do not run two scripts that write to the same folder at once. Later experiments read the outputs of earlier ones, so keep the upstream folders.

## Protocol notes

`CAPRAZ_TEKRAR_PROTOKOLU.md` and `FARKLI_PARAMETRELER_PROTOKOLU.md` are the plans for experiments 1 and 2, written in Turkish before the runs. They are left unchanged because the scripts hash them, so the regime protocol still ends with its pre-run status line. Scripts 3, 4 and 6 also write a Turkish draft `RAPOR.md` when rerun; that file is ignored by git, and `REPORT.md` is the edited version.
