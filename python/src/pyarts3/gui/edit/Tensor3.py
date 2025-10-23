"""Editor for Tensor3 (3D array) values."""

import numpy as np
from ..common import edit_ndarraylike

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit a Tensor3 (3D array) value.
    
    Parameters
    ----------
    value : Tensor3 or 3D array-like
        The 3D tensor to edit
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
