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
from . import ArrayOf
from . import Options
from . import UnifiedPropertyEditor


def get_editable_types():
    """
    Get a set of all type names that can be edited.
    
    With the UnifiedPropertyEditor, virtually all ARTS types can be edited:
    - Terminal types (Numeric, Index, String, etc.) use specialized editors
    - ArrayOf* types use the ArrayOf editor
    - Option enums use the Options editor
    - Types with properties use the UnifiedPropertyEditor
    - Array-like, GriddedField, and map-like types use specialized editors
    
    Returns all types from builtin_groups() since the unified editor can handle
    anything - it will show properties for complex types and route terminals
    to appropriate editors.
    
    Returns
    -------
    set of str
        Set of all ARTS type names (typically 348 types from builtin_groups)
    """
    editable = set()
    
    # Use builtin_groups() as the authoritative source of ARTS types
    # The unified editor can handle all of these
    try:
        from pyarts3.utils import builtin_groups
        builtin_types = builtin_groups()
        
        for type_class in builtin_types:
            type_name = type_class.__name__
            editable.add(type_name)
    except Exception:
        pass
    
    return editable


def can_edit(type_name):
    """
    Check if a type can be edited.
    
    Parameters
    ----------
    type_name : str
        The type name to check
        
    Returns
    -------
    bool
        True if the type can be edited
    """
    return type_name in get_editable_types()


def edit(value, parent=None):
    """
    Generic edit function that dispatches to the appropriate editor.

    Uses the unified property-based editor for ARTS types with properties,
    and delegates to specialized editors for terminal types.

    This function automatically determines the type of the input value and:
    1. For terminal types (Numeric, Index, String, ArrayOf*, Options, arrays, 
       Python built-ins), uses specialized editors
    2. For complex types with properties, uses the unified property editor
       with breadcrumb navigation

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

    >>> # Complex types use unified property editor
    >>> sun = pyarts.arts.Sun()
    >>> new_sun = edit.edit(sun)  # Opens breadcrumb property interface
    """
    # Special-case: Workspace must use its dedicated editor (never unified)
    type_name = type(value).__name__
    if type_name == 'Workspace':
        from . import Workspace as WorkspaceEditor
        return WorkspaceEditor.edit(value, parent=parent)

    # Use unified editor's logic for all routing
    from .UnifiedPropertyEditor import edit as unified_edit
    return unified_edit(value, parent=parent)
