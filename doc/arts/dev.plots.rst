.. _dev.plots:

Plot Modules for GUI Integration
==================================

Overview
--------

The ARTS Python interface includes a comprehensive set of plot modules for
visualizing workspace group data types. These modules are located in
``python/src/pyarts3/plots/`` and provide a consistent interface for creating
matplotlib visualizations.  It is non-exhaustive but covers several commonly
used workspace groups.

The intent is that we might use these plot modules in GUI applications
that need to
visualize ARTS data.  Each plot module is designed to handle a specific
workspace group type, providing tailored plotting functionality.

Design Pattern
--------------

All plot modules follow a consistent pattern with the ``plot()`` function:

.. code-block:: python

    def plot(
        data: pyarts.arts.<WorkspaceGroup>,
        *,
        fig: matplotlib.figure.Figure | None = None,
        ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
        # ... additional parameters
        **kwargs
    ) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
        """Plot the workspace group data.

        Parameters
        ----------
        data : ~pyarts3.arts.<WorkspaceGroup>
            The data to plot
        fig : ~matplotlib.figure.Figure, optional
            The matplotlib figure to draw on. Defaults to None for new figure.
        ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
            The matplotlib axes to draw on. Defaults to None for new axes.
        **kwargs
            Additional keyword arguments passed to underlying matplotlib functions

        Returns
        -------
        fig :
            As input if input.  Otherwise the created Figure.
        ax :
            As input if input.  Otherwise the created Axes.
        """
        # Implementation
        return fig, ax

**Critical:** All plot functions must accept ``**kwargs`` and pass them to the
underlying matplotlib plotting functions. This allows GUI frameworks to pass
arbitrary styling parameters without knowing the specific parameters each plot
type uses.

Usage Examples
--------------

Basic Usage
^^^^^^^^^^^

.. code-block:: python

    import pyarts3 as pyarts
    import pyarts3.plots as plots

    # Create or load data
    vector_data = pyarts.arts.Vector([1, 2, 3, 4, 5])

    # Create plot
    fig, ax = plots.Vector.plot(
        vector_data,
        xlabel="Index",
        ylabel="Value",
        title="My Vector Plot"
    )

    # Display or save
    import matplotlib.pyplot as plt
    plt.show()

Adding New Plot Modules
^^^^^^^^^^^^^^^^^^^^^^^

To add a new plot module for a workspace group:

1. Create ``python/src/pyarts3/plots/NewType.py`` following the standard
   pattern.
2. Implement the ``plot()`` function with appropriate visualization.
3. Add Sphinx documentation with ``.. plot::`` directive examples if possible.
4. Import the module in ``python/src/pyarts3/plots/__init__.py``

Required imports:

.. code-block:: python

    import pyarts3 as pyarts
    import matplotlib.pyplot as plt

    __all__ = ['plot']

Future Enhancements
-------------------

Please ensure to add new types as you implement them if you
deem them plot-worthy.
It helps ensure that we have a comprehensive set of
visualization tools for ARTS data,
and visualization is a key aspect of helping users
understand what we are doing.
