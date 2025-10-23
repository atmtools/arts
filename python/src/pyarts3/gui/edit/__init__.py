"""Editor modules for ARTS workspace group types.

This module provides editing capabilities for ARTS workspace types.
Each submodule is named exactly as the workspace group type it can edit.

Each submodule provides an edit() function that takes the data value as its
first argument and optional parent widget. The edit() functions return the
edited value if accepted, or None if cancelled.

Example usage
-------------
>>> import pyarts3.gui.edit as edit
>>> vec = pyarts.arts.Vector([1, 2, 3])
>>> new_vec = edit.edit(vec)

Or directly:
>>> from pyarts3.gui.edit import Vector
>>> new_vec = Vector.edit(vec)
"""

# Import available editor modules
from . import Generic
from . import Index
from . import Matrix
from . import Numeric
from . import Stokvec
from . import String
from . import Vector
from . import Vector2
from . import Vector3
from . import Workspace


def edit(value, parent=None):
    """
    Generic edit function that dispatches to the appropriate editor module.

    This function automatically determines the type of the input value and calls
    the corresponding editor module's edit() function.

    Parameters
    ----------
    value : ARTS workspace type
        Any ARTS data type that has a corresponding editor module.
    parent : QWidget, optional
        Parent widget for the dialog. Defaults to None.

    Returns
    -------
    edited_value or None
        The edited value if accepted, None if cancelled or read-only.

    Examples
    --------
    >>> import pyarts3 as pyarts
    >>> import pyarts3.gui.edit as edit
    >>> vec = pyarts.arts.Vector([1, 2, 3, 4])
    >>> new_vec = edit.edit(vec)

    >>> num = pyarts.arts.Numeric(3.14)
    >>> new_num = edit.edit(num)
    """
    import sys

    # Get the type name of the value
    type_name = type(value).__name__

    # Get the current module (pyarts3.gui.edit)
    current_module = sys.modules[__name__]

    # Check if we have a submodule with the same name as the type
    if hasattr(current_module, type_name):
        edit_module = getattr(current_module, type_name)

        # Check if the submodule has an edit function
        if hasattr(edit_module, 'edit'):
            return edit_module.edit(value, parent=parent)
        else:
            # Shouldn't happen, but fallback to Generic
            return Generic.edit(value, parent=parent)
    else:
        # No specific editor found, use generic viewer
        return Generic.edit(value, parent=parent)
