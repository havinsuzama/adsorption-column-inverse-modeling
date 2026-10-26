# Testing and environment

The archived experiments were run on Windows 10 with Python 3.14.6, NumPy 2.5.1, SciPy 1.18.0 and Matplotlib 3.11.1, and PyTorch 2.13.0 (CPU) for the PINN. The GitHub Actions workflow (`.github/workflows/checks.yml`) runs on Python 3.11.

## Checks

```bash
python -m unittest discover -s tests -v
```

- Mass balance and concentration bounds of the reference solver.
- Analytical parameter sensitivities against central finite differences.

```bash
python run_demo.py
python run_demo.py --mode results
python run_demo.py --mode inverse
python experimental_pinn/balanced_pinn.py --mode smoke --no-plot
```

The PINN smoke run only checks that the training loop runs; it does not evaluate the physics. None of these checks validate the model against experiments.

## Line endings

`.gitattributes` stores the files under `experiments/` and `experimental_pinn/` byte for byte. The experiment scripts check the SHA-256 hashes of their sources before resuming a run, so automatic line-ending conversion would make them refuse to resume.
