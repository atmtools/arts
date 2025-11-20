"""  Plotting ARTS data.

This module provides functions related to plotting ARTS data
based only on the type of the data.
Each submodule is named exactly as the data type it can plot, which is also
the name of a builtin group.  It is supposed to be a helpful way to get
a quick overview of the data in a given ARTS data structure.

Each submodule provides a plot()-method.
The generic plot() function in this module dispatches to the appropriate
submodule based on the type of the input data.  Only exact type matches are
considered.  For unavailable types, a TypeError is raised.
Each submodule's plot() function accepts the same 3 parameters:

- ``data``: The ARTS data to plot.
- ``fig``: An optional matplotlib Figure to plot on.
- ``ax``: An optional matplotlib Axes (or list/array of Axes) to plot on.

Additional keyword arguments are forwarded to the specific plot() function,
which if they are not recognized directly by the specific plot() function are
forwarded yet again to the underlying matplotlib plotting functions.

.. tip::
    Below you will find examples using mock-data to illustrate the usage of the
    specific plot-functions.  In practice, you would use ARTS data created by
    running ARTS simulations.
"""

import numpy as _numpy
import matplotlib as _matplotlib
from . import AbsorptionBands
from . import ArrayOfCIARecord
from . import ArrayOfPropagationPathPoint
from . import ArrayOfSensorObsel
from . import AscendingGrid
from . import AtmField
from . import AziGrid
from . import CIARecord
from . import DisortFlux
from . import DisortRadiance
from . import GeodeticField2
from . import GriddedField2
from . import LatGrid
from . import LonGrid
from . import MagneticAngles
from . import Matrix
from . import MuelmatVector
from . import PropmatMatrix
from . import PropmatVector
from . import SortedGriddedField1
from . import SortedGriddedField2
from . import StokvecMatrix
from . import StokvecVector
from . import SubsurfaceField
from . import Sun
from . import SurfaceField
from . import Vector
from . import ZenGrid


def plot(data: object,
         *,
         fig: _matplotlib.figure.Figure | None = None,
         ax: _matplotlib.axes.Axes | list[_matplotlib.axes.Axes] | _numpy.ndarray[_matplotlib.axes.Axes] | None = None,
         **kwargs) -> tuple[_matplotlib.figure.Figure, _matplotlib.axes.Axes | list[_matplotlib.axes.Axes] | _numpy.ndarray[_matplotlib.axes.Axes]]:
    """Generic plotting function that dispatches to the appropriate plot module.

    This function automatically determines the type of the input data and calls
    the corresponding plot module's plot() function. All keyword arguments are
    forwarded to the specific plotting function.

    For convenience, this function can be used directly from the main pyarts3 module:

    Examples
    --------
    >>> import pyarts3 as pyarts
    >>> vec = pyarts.arts.Vector([1, 2, 3, 4])
    >>> fig, ax = pyarts.plot(vec)

    >>> mat = pyarts.arts.Matrix([[1, 2], [3, 4]])
    >>> fig, ax = pyarts.plot(mat, cmap='viridis')

    Parameters
    ----------
    data : ARTS builtin group
        Any ARTS data type that has a corresponding plot module.
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.

    Raises
    ------
    TypeError
        If no plot module exists for the given data type.
    Exception
        As forwarded by the specific plot module.
    """

    import sys

    # Get the type name of the data
    type_name = type(data).__name__

    # Get the current module (pyarts3.plots)
    current_module = sys.modules[__name__]

    # Check if we have a submodule with the same name as the type
    if hasattr(current_module, type_name):
        return getattr(current_module, type_name).plot(data, fig=fig, ax=ax, **kwargs)
    else:
        available_modules = {', '.join([name for name in dir(current_module) if not name.startswith('_') and name != 'plot' and name != 'sys' and name != 'common'])}
        import pyarts3 as pyarts
        if type_name in [x.__name__ for x in pyarts.utils.builtin_groups()]:
            raise TypeError(
            f"No plot module found for type '{type_name}'.\n"
            "If you create a plotting routine for this type, please consider contributing it to pyarts3!\n"
            f"Available plot modules exist for: {available_modules}"
            )

        raise TypeError(
            f"No plot module found for type '{type_name}'.\n"
            f"Available plot modules exist for: {available_modules}"
        )
