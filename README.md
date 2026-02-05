# AE 246 Homework #1 (Python)

Clean, reusable implementation of strip-theory wing coefficients and wing-only pitching moment for AE 246 HW1.

## What’s included

- **Q1**: Compute `CL0`, `CLalpha`, `CL`, `CD`, and induced/downwash angle `epsilon = K*CL`
- **Q2**: Wing-only pitching moment about the CG: `Cm0`, `Cma`, `Cm`, and trim `alpha_trim`
- **Q3**: Compare linear vs cubic-nonlinear `Cm(alpha)` and plot for `0° ≤ alpha ≤ 10°`

## Repo structure

```
ae246_hw1/
  data/
    wing_geometry.m
  scripts/
    run_q1_q2.py
    run_q3.py
  src/
    ae246_hw1/
      __init__.py
      geometry.py
      aero.py
      moments.py
```

## Setup

```bash
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate

pip install -r requirements.txt
```

## Run

From the repo root:

```bash
python -m pip install -e .
python scripts/run_q1_q2.py --enforce-area
python scripts/run_q3.py --enforce-area
```

## Notes

- `--enforce-area` scales the tabulated chord so that the integrated area matches the given `S=1950 ft^2`.
- Angles inside formulas are **radians**; scripts accept degrees where appropriate.
