"""  Plotting ARTS data

This module provides functions related to plotting ARTS data
based only on the type of the data.
Each submodule is named exactly as the data type it can plot, which is also
the name of a builtin group.

Each submodule provides a plot()
function that takes the data type as its first argument and optional
arguments to control the plotting.  The plot() functions return the
figure, a list of subplots, and nothing more.
"""

from . import AbsorptionBands
from . import AtmField
from . import ArrayOfPropagationPathPoint
from . import ArrayOfSensorObsel
from . import MagneticAngles
from . import MuelmatVector
from . import SubsurfaceField
