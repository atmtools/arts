Plot GUI architecture
======================

Overview
--------

The GUI in ``python/src/pyarts3/gui/qt_plot_gui.py`` provides an interactive
front-end for plotting ARTS workspace data using the functions in
``pyarts3.plots`` and editing values using ``pyarts3.gui.edit``.

It is designed for experimentation: developers can wire in a plot function,
provide a callable that computes plotting keyword arguments, and interactively
adjust simulation settings, results, and additional options.

Main components
---------------

- ``PlotGui``: the main Qt widget hosting the toolbar, canvas, and three lists:
  - Simulation Settings
  - Results
  - Additional Options
- ``start_gui``: helper to create a QApplication (if needed) and
  show ``PlotGui``.
- Editors: double-click on any list item to edit its value via
  ``pyarts3.gui.edit.edit(value, parent)``.
- Workspace dialog: a dedicated editor lists workspace variables, supports
  adding unset variables, and dispatches to per-type editors.

Flow
----

1. Provide a plot function, e.g., ``pyarts3.plots.Vector.plot``.
2. Provide a callable that returns default keyword arguments for the plot
   (Simulation Settings → plot kwargs → Results are computed and displayed).
3. Edit values by double-clicking items in the lists; the plot re-renders after
   accepted edits.
4. Use Additional Options to pass arbitrary keyword arguments to the plot
   function.

Editors integration
-------------------

- The dispatcher in ``pyarts3.gui.edit`` resolves a per-type editor module by
  matching the value's type name. If missing, the ``Generic`` editor falls back
  to a type-preserving eval.
- Editors must preserve ARTS types (e.g., ``Vector``, ``Matrix``). Reconstruct
  using ``type(value)(payload)`` on save.

Extending the GUI
-----------------

- Add new plot modules under ``python/src/pyarts3/plots/`` (see ``dev.plots``).
- Add new editors under ``python/src/pyarts3/gui/edit/`` (see ``dev.edit``).
- Wire custom post-callbacks or additional UI affordances in
  ``qt_plot_gui.py``.
