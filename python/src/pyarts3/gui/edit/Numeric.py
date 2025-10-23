"""Editor for Numeric (floating point) values."""

from pyarts3.gui.common import edit_numeric

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit a Numeric (floating point) value.
    
    Parameters
    ----------
    value : Numeric or float
        The numeric value to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    float or None
        The edited value if accepted, None if cancelled
    """
    return edit_numeric(value, parent)
