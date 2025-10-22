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
from . import ArrayOfPropagationPathPoint
from . import ArrayOfSensorObsel
from . import AscendingGrid
from . import AtmField
from . import AtmPoint
from . import AzimuthGrid
from . import DisortFlux
from . import DisortRadiance
from . import GeodeticField2
from . import GriddedField2
from . import LatGrid
from . import LonGrid
from . import MagneticAngles
from . import Matrix
from . import MuelmatVector
from . import PropagationPathPoint
from . import PropmatMatrix
from . import PropmatVector
from . import SortedGriddedField1
from . import SortedGriddedField2
from . import SortedGriddedField3
from . import Stokvec
from . import StokvecMatrix
from . import StokvecVector
from . import SubsurfaceField
from . import Sun
from . import SurfaceField
from . import Tensor3
from . import Tensor4
from . import Time
from . import Vector
from . import ZenithGrid
