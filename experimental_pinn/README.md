# Experimental PINN (early-time window)

`balanced_pinn.py` fine-tunes a PINN for the start of the run (normalized time ≤ 0.01), beginning from `start_model.pt`. Collocation points are split evenly over four log-time bins so that the earliest times are not under-sampled. This is a diagnostic, not a working surrogate for the full breakthrough curve.

Run from the repository root:

```bash
python -m pip install -r requirements-pinn.txt
python experimental_pinn/balanced_pinn.py --mode smoke --no-plot
```

- `start_model.pt` is a checkpoint from earlier training runs in this project, loaded with `weights_only=True`.
- `reference_early.npz` holds finite-volume fields (`ref_c`, `ref_s`) used only by `--validate`, never for training.
- Every run starts from the same checkpoint, so repeated runs are not independent initializations.
- Output goes to `experimental_pinn/runs/`, which is not tracked by git.

The `--validate` gate has not been passed yet. `quick` and `standard` modes take noticeably longer than `smoke`.
