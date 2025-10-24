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


def plot(data, *, fig=None, ax=None, **kwargs):
    """Generic plotting function that dispatches to the appropriate plot module.

    This function automatically determines the type of the input data and calls
    the corresponding plot module's plot() function. All keyword arguments are
    forwarded to the specific plotting function.

    Parameters
    ----------
    data : ARTS workspace type
        Any ARTS data type that has a corresponding plot module.
    fig : matplotlib Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : matplotlib Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    **kwargs
        All keyword arguments are forwarded to the specific plot function.
        Additional plot-specific parameters vary by data type.

    Returns
    -------
    fig : matplotlib Figure
        The matplotlib figure object.
    ax : matplotlib Axes or list of Axes
        The matplotlib axes object(s).

    Raises
    ------
    TypeError
        If no plot module exists for the given data type.

    Examples
    --------
    >>> import pyarts3 as pyarts
    >>> vec = pyarts.arts.Vector([1, 2, 3, 4])
    >>> fig, ax = pyarts.plots.plot(vec)

    >>> mat = pyarts.arts.Matrix([[1, 2], [3, 4]])
    >>> fig, ax = pyarts.plots.plot(mat, cmap='viridis')
    """

    import sys


    # Get the type name of the data
    type_name = type(data).__name__

    # Get the current module (pyarts3.plots)
    current_module = sys.modules[__name__]

    # Check if we have a submodule with the same name as the type
    if hasattr(current_module, type_name):
        plot_module = getattr(current_module, type_name)

        # Check if the submodule has a plot function
        if hasattr(plot_module, 'plot'):
            return plot_module.plot(data, fig=fig, ax=ax, **kwargs)
        else:
            raise TypeError(
                f"Plot module '{type_name}' exists but has no plot() function")
    else:
        raise TypeError(
            f"No plot module found for type '{type_name}'. "
            f"Available plot modules: {', '.join([name for name in dir(current_module) if not name.startswith('_') and name != 'plot' and name != 'sys'])}"
        )
