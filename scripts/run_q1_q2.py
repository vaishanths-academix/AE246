from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ae246_hw1.geometry import load_wing_geometry_m
from ae246_hw1.aero import WingAeroInputs, WingStripTheory
from ae246_hw1.moments import MomentInputs, wing_cm_linear


def main() -> None:
    p = argparse.ArgumentParser(description="AE246 HW1 - Q1/Q2 quick runner")
    p.add_argument("--geom", type=Path, default=Path("data/wing_geometry.m"), help="Path to wing_geometry.m")
    p.add_argument("--alpha-deg", type=float, default=1.9017, help="Angle of attack (deg)")
    p.add_argument("--enforce-area", action="store_true", help="Scale chord so geometry area matches S")

    args = p.parse_args()

    # Inputs from homework statement
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

    res = wing.analyze(alpha_deg=args.alpha_deg)

    print("\n=== Q1: Aerodynamics ===")
    print(f"alpha = {args.alpha_deg:.4f} deg ({np.deg2rad(args.alpha_deg):.6f} rad)")
    print(f"CL0        = {res.CL0:.8f}")
    print(f"CLalpha    = {res.CLalpha:.8f} 1/rad")
    print(f"CL         = {res.CL:.8f}")
    print(f"CD         = {res.CD:.8f}")
    print(f"CDi        = {res.CDi:.8f}")
    print(f"epsilon_i  = {res.epsilon_ind_deg:.6f} deg  (epsilon_i = K*CL)")
    print(f"b          = {res.b:.4f} ft")
    print(f"AR         = {res.AR:.6f}")
    print(f"S(raw)     = {res.S_raw_from_data:.4f} ft^2")
    print(f"chordScale = {res.chord_scale:.8f}")

    # Q2 moment about CG (wing-only, Cmac = 0)
    m = MomentInputs(
        cbar=inputs.cbar,
        xcg_from_le_mac=3.8276,
        xac_from_le_mac=0.25 * inputs.cbar,
        Cmac=0.0,
    )
    alpha_rad = float(np.deg2rad(args.alpha_deg))
    mres = wing_cm_linear(res.CL0, res.CLalpha, alpha_rad, m)

    print("\n=== Q2: Wing pitching moment about CG (linear) ===")
    print(f"delta = (xcg-xac)/cbar = {m.delta:.6f}")
    print(f"Cm0   = {mres.Cm0:.10f}")
    print(f"Cma   = {mres.Cma:.10f} 1/rad")
    print(f"Cm    = {mres.Cm:.10f}")
    print(f"alpha_trim = {mres.alpha_trim_deg:.6f} deg ({mres.alpha_trim_rad:.6f} rad)")


if __name__ == "__main__":
    main()
