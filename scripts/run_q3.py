from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from ae246_hw1.geometry import load_wing_geometry_m
from ae246_hw1.aero import WingAeroInputs, WingStripTheory
from ae246_hw1.moments import MomentInputs, wing_cm_linear, cm_nonlinear_cubic_series, eval_cubic


def main() -> None:
    p = argparse.ArgumentParser(description="AE246 HW1 - Q3(b,c): Cm(10deg) and plot")
    p.add_argument("--geom", type=Path, default=Path("data/wing_geometry.m"), help="Path to wing_geometry.m")
    p.add_argument("--enforce-area", action="store_true", help="Scale chord so geometry area matches S")
    p.add_argument("--alpha-max-deg", type=float, default=10.0, help="Max alpha for plot (deg)")
    args = p.parse_args()

    inputs = WingAeroInputs(
        S=1950.0,
        cbar=16.6417,
        cl0_2d=0.1042,
        cla_2d=7.2856,
        cd0=0.0208,
        K=0.0474,
    )

    geom = load_wing_geometry_m(args.geom, incidence_units="deg")
    wing = WingStripTheory(geom, inputs, enforce_area=args.enforce_area)

    # We need CL0 and CLalpha from the wing model (independent of alpha value)
    CL0, CLalpha, _ = wing.CL_components(alpha_rad=0.0)

    m = MomentInputs(
        cbar=inputs.cbar,
        xcg_from_le_mac=3.8276,
        xac_from_le_mac=0.25 * inputs.cbar,
        Cmac=0.0,
    )

    # Linear Cm model coefficients
    lin = wing_cm_linear(CL0, CLalpha, alpha_rad=0.0, m=m)

    # Nonlinear cubic coefficients
    coeffs = cm_nonlinear_cubic_series(CL0, CLalpha, inputs.cd0, inputs.K, m)

    # (b) Cm at alpha=10 deg
    a10 = float(np.deg2rad(10.0))
    Cm_lin_10 = lin.Cm0 + lin.Cma * a10
    Cm_nl_10 = eval_cubic(a10, coeffs)

    print("\n=== Q3(b): Cm at alpha = 10 deg ===")
    print(f"Cm_linear(10 deg)   = {Cm_lin_10:.10f}")
    print(f"Cm_nonlinear(10 deg)= {Cm_nl_10:.10f}")
    print("\nNonlinear cubic polynomial (alpha in rad):")
    a0,a1,a2,a3 = coeffs
    print(f"Cm(a) = {a0:.6e} + {a1:.6e} a + {a2:.6e} a^2 + {a3:.6e} a^3")

    # (c) plot 0..alpha_max
    alpha_deg = np.linspace(0.0, args.alpha_max_deg, 401)
    alpha_rad = np.deg2rad(alpha_deg)

    Cm_lin = lin.Cm0 + lin.Cma * alpha_rad
    Cm_nl = eval_cubic(alpha_rad, coeffs) if hasattr(alpha_rad, "__len__") else eval_cubic(float(alpha_rad), coeffs)

    # vectorize eval_cubic for arrays
    Cm_nl = a0 + a1*alpha_rad + a2*alpha_rad**2 + a3*alpha_rad**3

    plt.figure()
    plt.plot(alpha_deg, Cm_lin, label="Linear $C_m$")
    plt.plot(alpha_deg, Cm_nl, label="Nonlinear cubic $C_m$")
    plt.xlabel(r"$\alpha$ (deg)")
    plt.ylabel(r"$C_m$")
    plt.title(r"$C_m(\alpha)$, Wing Alone")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
