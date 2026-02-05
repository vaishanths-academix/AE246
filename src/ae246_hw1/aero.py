from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .geometry import WingGeometry


@dataclass(frozen=True)
class WingAeroInputs:
    """Aerodynamic and planform inputs."""
    S: float                 # wing area used for normalization [ft^2]
    cbar: float              # mean aerodynamic chord [ft]
    cl0_2d: float            # section cl0
    cla_2d: float            # section cl_alpha [1/rad]
    cd0: float               # zero-lift drag coefficient
    K: float                 # induced drag factor


@dataclass(frozen=True)
class WingAeroResults:
    """Computed aerodynamic coefficients."""
    CL0: float
    CLalpha: float
    CL: float
    CD: float
    CDi: float
    epsilon_ind_rad: float
    epsilon_ind_deg: float
    AR: float
    b: float
    S_raw_from_data: float
    chord_scale: float


class WingStripTheory:
    """
    Wing strip-theory integrator for CL using tabulated geometry.

    Assumptions:
      - Section cl0 and cl_alpha are constant along the span.
      - Effective section angle = alpha + i(y), i in radians.
      - Symmetric wing: integrate over semi-span and multiply by 2.
    """

    def __init__(self, geom: WingGeometry, inputs: WingAeroInputs, enforce_area: bool = True) -> None:
        self.geom = geom
        self.inputs = inputs
        self.enforce_area = enforce_area

    def _scaled_chord(self) -> tuple[np.ndarray, float, float]:
        """
        Returns chord scaled so that geometry-integrated area equals inputs.S (optional).
        Returns: (chord_used, S_raw, scale)
        """
        S_raw = self.geom.area()
        scale = 1.0
        chord = self.geom.chord.copy()
        if self.enforce_area:
            scale = self.inputs.S / S_raw
            chord *= scale
        return chord, S_raw, scale

    def CL_components(self, alpha_rad: float) -> tuple[float, float, float]:
        """
        Returns (CL0, CLalpha, CL) where:
          CL0 = CL at alpha=0 (still includes twist/incidence)
          CLalpha = dCL/dalpha (per rad)
          CL = CL(alpha_rad)
        """
        chord, _, _ = self._scaled_chord()
        y = self.geom.y
        i = self.geom.incidence_rad

        S = self.inputs.S
        cl0 = self.inputs.cl0_2d
        cla = self.inputs.cla_2d

        # CL0
        CL0 = (2.0 / S) * float(np.trapezoid((cl0 + cla * i) * chord, y))

        # CLalpha (should equal cla if area normalization is perfect)
        CLalpha = (2.0 / S) * float(np.trapezoid((cla) * chord, y))

        CL = CL0 + CLalpha * float(alpha_rad)
        return CL0, CLalpha, CL

    def drag_and_downwash(self, CL: float) -> tuple[float, float, float, float]:
        """
        Returns (CD, CDi, epsilon_ind_rad, epsilon_ind_deg)
        Uses polar CD = CD0 + K CL^2 and induced angle epsilon = K CL (rad).
        """
        CDi = self.inputs.K * CL**2
        CD = self.inputs.cd0 + CDi
        eps_rad = self.inputs.K * CL
        eps_deg = float(np.degrees(eps_rad))
        return float(CD), float(CDi), float(eps_rad), eps_deg

    def analyze(self, alpha_deg: float) -> WingAeroResults:
        """Convenience end-to-end analysis at a given alpha in degrees."""
        alpha_rad = float(np.deg2rad(alpha_deg))

        chord, S_raw, scale = self._scaled_chord()
        b = self.geom.span
        AR = b**2 / self.inputs.S

        CL0, CLalpha, CL = self.CL_components(alpha_rad)
        CD, CDi, eps_rad, eps_deg = self.drag_and_downwash(CL)

        return WingAeroResults(
            CL0=CL0,
            CLalpha=CLalpha,
            CL=CL,
            CD=CD,
            CDi=CDi,
            epsilon_ind_rad=eps_rad,
            epsilon_ind_deg=eps_deg,
            AR=float(AR),
            b=float(b),
            S_raw_from_data=float(S_raw),
            chord_scale=float(scale),
        )
