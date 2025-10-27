"""GUI subpackage for pyarts3.

Exports:
- start_gui: Convenience function to start the GUI
- PlotGui: Main plotting GUI widget
- edit: Module containing workspace group editors
"""
from .qt_plot_gui import PlotGui, start_gui
from . import edit

__all__ = ["start_gui", "PlotGui", "edit"]

print ("WARNING: Using this GUI module is highly experimental and may lead to crashes.")
