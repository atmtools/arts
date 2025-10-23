"""Editor for Matrix (2D array) values."""

import numpy as np
from ..common import edit_ndarraylike

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit a Matrix (2D array) value.
    
    Parameters
    ----------
    value : Matrix or 2D array-like
        The matrix value to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    numpy.ndarray or None
        The edited value if accepted, None if cancelled
    """
    result = edit_ndarraylike(value, parent)
    if result is None:
        return None
    
    # Preserve original ARTS type if possible
    try:
        # Convert to nested list for ARTS constructor
        result_list = result.tolist()
        return type(value)(result_list)
    except Exception:
        # Fallback to numpy array
        return result
