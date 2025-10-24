"""Editor for String values."""

from .widgets import edit_string

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit a String value.
    
    Parameters
    ----------
    value : String or str
        The string value to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    str or None
        The edited value if accepted, None if cancelled
    """
    return edit_string(value, parent)
