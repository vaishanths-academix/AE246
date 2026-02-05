from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class MomentInputs:
    """Moment geometry inputs about CG."""
    cbar: float        # MAC [ft]
    xcg_from_le_mac: float  # CG measured aft of LE_MAC [ft]
    xac_from_le_mac: float  # AC measured aft of LE_MAC [ft]
    Cmac: float = 0.0  # wing moment about AC (given as 0 in HW)

    @property
    def delta(self) -> float:
        """
        Nondimensional arm (xcg - xac)/cbar.
        Positive means CG aft of AC.
        """
        return (self.xcg_from_le_mac - self.xac_from_le_mac) / self.cbar


@dataclass(frozen=True)
class MomentResults:
    Cm0: float
    Cma: float
    Cm: float
    alpha_trim_rad: float
    alpha_trim_deg: float


def wing_cm_linear(CL0: float, CLalpha: float, alpha_rad: float, m: MomentInputs) -> MomentResults:
    """
    Linear wing pitching moment about CG:
      Cm = Cmac + (xcg-xac)/cbar * CL
    where CL = CL0 + CLalpha * alpha.
    """
    CL = CL0 + CLalpha * alpha_rad
    Cm0 = m.Cmac + m.delta * CL0
    Cma = m.delta * CLalpha
    Cm = Cm0 + Cma * alpha_rad

    # Trim: Cm=0 => alpha = -Cm0/Cma (if Cma != 0)
    if abs(Cma) < 1e-12:
        alpha_trim = float("nan")
    else:
        alpha_trim = -Cm0 / Cma

    return MomentResults(
        Cm0=float(Cm0),
        Cma=float(Cma),
        Cm=float(Cm),
        alpha_trim_rad=float(alpha_trim),
        alpha_trim_deg=float(np.degrees(alpha_trim)),
    )


def cm_nonlinear_cubic_series(
    CL0: float,
    CLalpha: float,
    CD0: float,
    K: float,
    m: MomentInputs,
) -> tuple[float, float, float, float]:
    """
    Build a cubic polynomial approximation for Cm(alpha) (alpha in rad) using:
      Cm = delta * [ CL cos(a) + CD sin(a) ] + Cmac
    with:
      CL = CL0 + CLalpha a
      CD = CD0 + K CL^2
    and series:
      cos(a) ≈ 1 - a^2/2
      sin(a) ≈ a - a^3/6
    Returns coefficients (a0,a1,a2,a3) for:
      Cm(a) ≈ a0 + a1 a + a2 a^2 + a3 a^3
    """
    delta = m.delta

    D0 = CD0 + K * CL0**2
    D1 = 2.0 * K * CL0 * CLalpha
    D2 = K * CLalpha**2

    a0 = m.Cmac + delta * (CL0)
    a1 = delta * (CLalpha + D0)
    a2 = delta * (D1 - 0.5 * CL0)
    a3 = delta * (D2 - 0.5 * CLalpha - (1.0 / 6.0) * D0)
    return float(a0), float(a1), float(a2), float(a3)


def eval_cubic(a: float, coeffs: tuple[float, float, float, float]) -> float:
    a0, a1, a2, a3 = coeffs
    return a0 + a1 * a + a2 * a**2 + a3 * a**3
