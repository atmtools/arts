"""Editor for Index (integer) values."""

from pyarts3.gui.common import edit_index

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit an Index (integer) value.
    
    Parameters
    ----------
    value : Index or int
        The integer value to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    int or None
        The edited value if accepted, None if cancelled
    """
    return edit_index(value, parent)
