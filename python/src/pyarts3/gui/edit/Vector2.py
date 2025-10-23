"""Editor for Vector2 (2D vector) values."""

import numpy as np
from ..common import edit_ndarraylike

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit a Vector2 (2D vector) value.
    
    Parameters
    ----------
    value : Vector2 or 2-element array
        The 2D vector to edit (X, Y components)
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
        # Convert to list for ARTS constructor
        result_list = result.tolist()
        return type(value)(result_list)
    except Exception:
        return result
