"""
AE 246 HW#1 — Q1, Q2, Q3(b,c)

@Author: Vaishanth Srinivasan

Assumes wing_geometry.m contains:
    data = [ y   chord   incidence_deg ];

Internal units:
- angles in radians
- lengths in ft
- outputs in both rad/deg where helpful
"""

from __future__ import annotations

import re
import math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


# =========================
# Constants
# =========================
S = 1950.0           # ft^2
cbar = 16.6417       # ft (MAC)
dihedral = 28.4286 * math.pi/180

cl0_2d = 0.1042
cla_2d = 7.2856       # 1/rad

alpha_q1_deg = 1.9017  # deg (for Q1a evaluation)
CD0 = 0.0208
K = 0.0474

# Q2 geometry: along fuselage centerline measured from LE of MAC
xcg_from_LEmac = 3.8276  # ft (CG behind LE of MAC)
xac_from_LEmac = 0.25 * cbar  # ft (AC at quarter-chord of MAC)


def parse_wing_geometry_m(path: str | Path = "wing_geometry.m") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse wing_geometry.m expecting:
        data = [ y   chord   incidence_deg ];
    Returns:
        y (ft), c (ft), i (rad)
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Could not find {path.resolve()}. Put wing_geometry.m in the same folder as this script."
        )

    text = path.read_text(errors="ignore")

    m = re.search(r"data\s*=\s*\[(.*?)\];", text, flags=re.S)
    if not m:
        raise ValueError("Could not find a 'data = [ ... ];' block in wing_geometry.m")

    block = m.group(1)
    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", block)
    vals = np.array([float(x) for x in nums], dtype=float)

    if vals.size % 3 != 0:
        raise ValueError(
            f"Expected 3 columns (y, chord, incidence). Found {vals.size} numbers (not divisible by 3)."
        )

    data = vals.reshape((-1, 3))
    y = data[:, 0]
    c = data[:, 1]
    i_deg = data[:, 2]
    i = np.deg2rad(i_deg)

    return y, c, i


def enforce_area(y: np.ndarray, c: np.ndarray, S_target: float) -> tuple[np.ndarray, float, float]:
    """
    Scales chord uniformly so that 2*∫c(y)dy matches S_target.
    Returns scaled chord, scale factor, and resulting S_check.
    """
    S_half_raw = np.trapezoid(c, y)
    S_raw = 2.0 * S_half_raw
    scale = S_target / S_raw
    c_scaled = c * scale
    S_check = 2.0 * np.trapezoid(c_scaled, y)
    return c_scaled, scale, S_check


# =========================
# Q1: Lift + drag + induced/downwash angle
# =========================
def solve_q1(y: np.ndarray, c: np.ndarray, i: np.ndarray) -> dict:
    """
    Computes CL0, CLalpha, CL at alpha_q1_deg; CD; AR; downwash epsilon.
    """
    alpha = np.deg2rad(alpha_q1_deg)

    # ---- Correct swept-wing CL0 ----
    # Camber term
    CL0_camber = cl0_2d * np.cos(dihedral) ** 2

    # Twist term
    twist_integral = np.trapezoid(i * c, y)
    CL0_twist = (2.0 * cla_2d * np.cos(dihedral) / S) * twist_integral

    CL0 = CL0_camber + CL0_twist

    # CLalpha = (2/S) ∫ cla * c(y) dy
    CLalpha = cla_2d * math.cos(dihedral)

    CL = CL0 + CLalpha * alpha

    # Drag polar
    CD = CD0 + K * CL**2
    CDi = K * CL**2

    # Geometry-derived AR
    b = 2.0 * float(np.max(y))
    AR = b**2 / S

    # Induced/downwash angle (lifting-line small-angle approximation):
    # epsilon ≈ CDi / CL = K*CL  (radians)
    epsilon_ind_rad = K * CL
    epsilon_ind_deg = np.degrees(epsilon_ind_rad)

    # Alternative "downwash" estimate some texts use: eps ≈ 2CL/(pi*AR)
    epsilon_alt_rad = 2.0 * CL / (np.pi * AR)
    epsilon_alt_deg = np.degrees(epsilon_alt_rad)

    return {
        "alpha_rad": alpha,
        "b": b,
        "AR": AR,
        "CL0": CL0,
        "CLalpha": CLalpha,
        "CL": CL,
        "CD": CD,
        "CDi": CDi,
        "eps_ind_rad": epsilon_ind_rad,
        "eps_ind_deg": epsilon_ind_deg,
        "eps_alt_rad": epsilon_alt_rad,
        "eps_alt_deg": epsilon_alt_deg,
    }


# =========================
# Q2: Wing pitching moment about CG (linear)
# =========================
def solve_q2(CL0: float, CLalpha: float, CL: float) -> dict:
    """
    Using Cmac = 0 and AC at 1/4 MAC:
      Cm = Cm_ac + CL * ((xcg/cbar) - (xac/cbar))
    Returns Cm0, Cma, Cm(at alpha_q1), alpha_trim.
    """
    Cm_ac = 0.0
    arm_nd = (xcg_from_LEmac / cbar) - (xac_from_LEmac / cbar)  # nondimensional moment arm

    Cm0_w = Cm_ac + CL0 * arm_nd
    Cma_w = CLalpha * arm_nd
    Cm_w = CL * arm_nd

    # Trim (wing-alone): set Cm = 0 => alpha_trim = -Cm0/Cma
    alpha_trim_rad = -Cm0_w / Cma_w
    alpha_trim_deg = np.degrees(alpha_trim_rad)

    return {
        "arm_nd": arm_nd,
        "Cm0_w": Cm0_w,
        "Cma_w": Cma_w,
        "Cm_w": Cm_w,
        "alpha_trim_rad": alpha_trim_rad,
        "alpha_trim_deg": alpha_trim_deg,
    }


# =========================
# Q3: Cm linear vs nonlinear polynomial (b,c)
# =========================
def q3_polynomial_coeffs(CL0: float, CLalpha: float) -> tuple[float, float, float, float]:
    """
    Builds a cubic polynomial for Cm(alpha):
      Cm(alpha) = a0 + a1*a + a2*a^2 + a3*a^3   (a in radians)

    Derived from:
      Cm = delta [ CL cos(a) + CD sin(a) ]
    with:
      CL = CL0 + CLalpha*a
      CD = CD0 + K*CL^2
      cos(a)≈1-a^2/2, sin(a)≈a-a^3/6
    and collecting terms up to a^3.

    delta = (xcg - xac)/cbar (nondimensional arm)
    """
    delta = (xcg_from_LEmac - xac_from_LEmac) / cbar

    D0 = CD0 + K * CL0**2
    D1 = 2.0 * K * CL0 * CLalpha
    D2 = K * CLalpha**2

    a0 = delta * (CL0)
    a1 = delta * (CLalpha + D0)
    a2 = delta * (D1 - 0.5 * CL0)
    a3 = delta * (D2 - 0.5 * CLalpha - (1.0 / 6.0) * D0)

    return a0, a1, a2, a3


def Cm_linear(alpha_rad: np.ndarray, Cm0: float, Cma: float) -> np.ndarray:
    return Cm0 + Cma * alpha_rad


def Cm_nonlinear(alpha_rad: np.ndarray, a0: float, a1: float, a2: float, a3: float) -> np.ndarray:
    return a0 + a1 * alpha_rad + a2 * alpha_rad**2 + a3 * alpha_rad**3


# =========================
# Main
# =========================
def main() -> None:
    # Load geometry
    y, c_raw, i = parse_wing_geometry_m("wing_geometry.m")

    # Enforce exact area
    c, chord_scale, S_check = enforce_area(y, c_raw, S_target=S)

    # Q1
    q1 = solve_q1(y, c, i)

    # Q2
    q2 = solve_q2(q1["CL0"], q1["CLalpha"], q1["CL"])

    # Q3 poly
    a0, a1, a2, a3 = q3_polynomial_coeffs(q1["CL0"], q1["CLalpha"])

    # Q3(b): alpha = 10 deg
    alpha10 = np.deg2rad(10.0)
    Cm_lin_10 = float(Cm_linear(alpha10, q2["Cm0_w"], q2["Cma_w"]))
    Cm_nl_10 = float(Cm_nonlinear(alpha10, a0, a1, a2, a3))

    # =========================
    # Print report
    # =========================
    print("\n" + "#" * 36)
    print("AE 246 HW#1 — Q1/Q2/Q3(b,c)")
    print("#" * 36 + "\n")

    print("== Geometry check ==")
    S_raw = 2.0 * np.trapezoid(c_raw, y)
    print(f"S raw from file:        {S_raw:.6f} ft^2")
    print(f"Chord scale applied:    {chord_scale:.10f}")
    print(f"S after scaling:        {S_check:.6f} ft^2")
    print(f"Span b:                 {q1['b']:.6f} ft")
    print(f"Aspect ratio AR:        {q1['AR']:.8f}\n")

    print("== Question 1 ==")
    print(f"CL0:                    {q1['CL0']:.8f}")
    print(f"CLalpha (1/rad):        {q1['CLalpha']:.8f}")
    print(f"CL @ {alpha_q1_deg:.4f} deg:        {q1['CL']:.8f}")
    print(f"CD:                     {q1['CD']:.8f}")
    print(f"CDi:                    {q1['CDi']:.8f}")
    print(f"epsilon_ind = K*CL:      {q1['eps_ind_rad']:.8e} rad = {q1['eps_ind_deg']:.6f} deg")
    print(f"epsilon_alt = 2CL/piAR:  {q1['eps_alt_rad']:.8e} rad = {q1['eps_alt_deg']:.6f} deg\n")

    print("== Question 2 ==")
    print(f"Moment arm (nd) (xcg/cbar - xac/cbar): {q2['arm_nd']:.8f}")
    print(f"Cm0_w:                  {q2['Cm0_w']:.10f}")
    print(f"Cma_w (1/rad):          {q2['Cma_w']:.10f}")
    print(f"Cm_w @ {alpha_q1_deg:.4f} deg:         {q2['Cm_w']:.10f}")
    print(f"alpha_trim:             {q2['alpha_trim_rad']:.10f} rad = {q2['alpha_trim_deg']:.6f} deg\n")

    print("== Question 3(b) ==")
    print("Cm at alpha = 10 deg:")
    print(f"  Linear:               {Cm_lin_10:.10f}")
    print(f"  Nonlinear (cubic):    {Cm_nl_10:.10f}\n")

    print("Nonlinear polynomial (alpha in radians):")
    print(f"  Cm(alpha) = {a0:.10e} + {a1:.10e}*a + {a2:.10e}*a^2 + {a3:.10e}*a^3\n")

    # =========================
    # Q3(c) plot
    # =========================
    alpha_deg = np.linspace(0.0, 10.0, 401)
    alpha_rad = np.deg2rad(alpha_deg)

    Cm_lin_curve = Cm_linear(alpha_rad, q2["Cm0_w"], q2["Cma_w"])
    Cm_nl_curve = Cm_nonlinear(alpha_rad, a0, a1, a2, a3)

    plt.figure()
    plt.plot(alpha_deg, Cm_lin_curve, label="Linear $C_m$")
    plt.plot(alpha_deg, Cm_nl_curve, label="Nonlinear cubic $C_m$")
    plt.xlabel(r"$\alpha$ (deg)")
    plt.ylabel(r"$C_m$")
    plt.title("Cm(alpha) Wing Alone, 0° ≤ α ≤ 10°")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
