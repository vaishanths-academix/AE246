from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np


@dataclass(frozen=True)
class WingGeometry:
    """Tabulated wing geometry (semi-span)."""
    y: np.ndarray          # spanwise stations, 0..b/2 [ft]
    chord: np.ndarray      # chord at stations [ft]
    incidence_rad: np.ndarray  # incidence/twist at stations [rad]

    @property
    def span(self) -> float:
        """Full wingspan b [ft], assuming y runs 0..b/2."""
        return 2.0 * float(np.max(self.y))

    def semi_area(self) -> float:
        """Semi-wing area from tabulated chord [ft^2]."""
        return float(np.trapezoid(self.chord, self.y))

    def area(self) -> float:
        """Full wing area from tabulated chord [ft^2]."""
        return 2.0 * self.semi_area()


def load_wing_geometry_m(path: str | Path, incidence_units: str = "deg") -> WingGeometry:
    """
    Load geometry from a MATLAB .m file containing:

        data = [ y   chord   incidence ];

    Parameters
    ----------
    path : str | Path
        Path to wing_geometry.m
    incidence_units : {"deg","rad"}
        Units of the incidence column in the file.

    Returns
    -------
    WingGeometry
    """
    path = Path(path)
    text = path.read_text(errors="ignore")

    m = re.search(r"data\s*=\s*\[(.*?)\];", text, flags=re.S)
    if not m:
        raise ValueError(f"Could not find a 'data = [ ... ];' block in {path}")

    block = m.group(1)
    nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", block)
    vals = np.array([float(x) for x in nums], dtype=float)

    if vals.size % 3 != 0:
        raise ValueError(
            f"Expected 3 columns (y, chord, incidence). Found {vals.size} numbers (not divisible by 3)."
        )

    data = vals.reshape((-1, 3))
    y = data[:, 0]
    chord = data[:, 1]
    inc = data[:, 2]

    units = incidence_units.lower().strip()
    if units == "deg":
        incidence_rad = np.deg2rad(inc)
    elif units == "rad":
        incidence_rad = inc
    else:
        raise ValueError("incidence_units must be 'deg' or 'rad'")

    return WingGeometry(y=y, chord=chord, incidence_rad=incidence_rad)
