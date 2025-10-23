"""Generic editor for array-like values (objects with __array__).

Uses the shared edit_ndarraylike() from common to provide a consistent UI
for 1D, 2D, and higher-dimensional arrays.
"""

import numpy as np
from ..common import edit_ndarraylike

__all__ = ["edit"]


def edit(value, parent=None):
    """
    Edit any array-like value (NumPy arrays or objects exposing __array__).

    Parameters
    ----------
    value : array-like
        The array value to edit (numpy arrays or classes exposing __array__).
    parent : QWidget, optional
        Parent widget for the dialog

    Returns
    -------
    array-like or None
        The edited value with original type preserved if accepted, None if cancelled
    """
    # edit_ndarraylike already handles type preservation, just return its result
    return edit_ndarraylike(value, parent)
