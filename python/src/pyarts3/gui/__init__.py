"""GUI subpackage for pyarts3.

Exports:
- start: Convenience function to start the GUI
"""
from .qt_plot_gui import PlotGui, start_gui

__all__ = ["start_gui"]
