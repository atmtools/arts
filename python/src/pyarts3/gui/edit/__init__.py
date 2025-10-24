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
from . import Time  # Treat Time as a basic type with its own simple editor
from . import UnifiedPropertyEditor


def get_editable_types():
    """
    Get a set of all type names that can be edited.
    
    Uses pyarts3.utils.builtin_groups() as the authoritative source of ARTS types,
    then determines which can be edited based on:
    - Explicit editor modules (Index, Numeric, String, ArrayOf, Options)
    - Option group enums (from arts.globals.option_groups())
    - Types with __array__ support (can use NDarray editor)
    - All ArrayOf* variants (can use ArrayOf editor)
    
    Returns
    -------
    set of str
        Set of all editable type names (typically 210+ types)
    """
    editable = set()
    
    # 1. Add explicit editor modules (excluding special ones)
    import sys
    current_module = sys.modules[__name__]
    for name in dir(current_module):
        if name.startswith('_'):
            continue
        if name in ('edit', 'get_editable_types', 'can_edit', 'Generic', 'Workspace', 'NDarray', 'ndarray'):
            continue
        attr = getattr(current_module, name)
        if hasattr(attr, 'edit'):
            editable.add(name)
    
    # 2. Add all option group enums
    try:
        import pyarts3.arts as arts
        option_groups = arts.globals.option_groups()
        editable.update(option_groups)
    except Exception:
        pass
    
    # 3. Use builtin_groups() as the authoritative source of ARTS types
    # This is the complete list of all ARTS builtin types (348 types)
    try:
        from pyarts3.utils import builtin_groups
        builtin_types = builtin_groups()
        
        for type_class in builtin_types:
            type_name = type_class.__name__
            
            # All ArrayOf* types can be edited via ArrayOf editor
            if type_name.startswith('ArrayOf'):
                editable.add(type_name)
                continue
            
            # Check if this type is a gridded field
            try:
                instance = type_class()
                if (hasattr(instance, '__array__') and 
                    hasattr(instance, 'grids') and 
                    hasattr(instance, 'gridnames') and 
                    hasattr(instance, 'dataname')):
                    editable.add(type_name)
                    continue
            except Exception:
                pass
            
            # Check if this type is map-like
            try:
                instance = type_class()
                if (hasattr(instance, 'keys') and 
                    hasattr(instance, 'items') and 
                    hasattr(instance, 'values') and
                    hasattr(instance, '__getitem__') and
                    hasattr(instance, '__setitem__')):
                    editable.add(type_name)
                    continue
            except Exception:
                pass
            
            # Check if this type has __array__ support (can use NDarray editor)
            if hasattr(type_class, '__array__'):
                editable.add(type_name)
                continue
            
            # Try to create an instance and check for __array__
            try:
                instance = type_class()
                if hasattr(instance, '__array__'):
                    editable.add(type_name)
            except Exception:
                pass  # Can't instantiate or no __array__, skip
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
       with tabbed navigation

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
    >>> new_sun = edit.edit(sun)  # Opens tabbed property interface
    """
    import sys
    
    # Special-case: Workspace must use its dedicated editor (never unified)
    type_name = type(value).__name__
    if type_name == 'Workspace':
        from . import Workspace as WorkspaceEditor
        return WorkspaceEditor.edit(value, parent=parent)

    # Try unified property editor first for complex ARTS types
    from .UnifiedPropertyEditor import edit as unified_edit, is_terminal_type, inspect_type
    
    # Get the type name of the value
    
    # For terminal types, use specialized editors
    if is_terminal_type(value):
        # Route to specialized editors
        current_module = sys.modules[__name__]
        
        # Check if we have a specific editor module
        if hasattr(current_module, type_name):
            edit_module = getattr(current_module, type_name)
            if hasattr(edit_module, 'edit'):
                return edit_module.edit(value, parent=parent)
        
        # Route ArrayOf* to the generic ArrayOf editor
        if type_name.startswith('ArrayOf'):
            return ArrayOf.edit(value, parent=parent)
        
        # Check if this is an option group enum
        import pyarts3.arts as arts
        option_groups = arts.globals.option_groups()
        if type_name in option_groups:
            return Options.edit(value, parent=parent)
        
        # Check if this is a gridded field
        is_griddedfield = (
            hasattr(value, '__array__') and
            hasattr(value, 'grids') and
            hasattr(value, 'gridnames') and
            hasattr(value, 'dataname')
        )
        
        if is_griddedfield:
            from ..common import edit_griddedfield
            return edit_griddedfield(value, parent=parent)
        
        # Check if this is a map-like object
        is_maplike = (
            hasattr(value, 'keys') and
            hasattr(value, 'items') and
            hasattr(value, 'values') and
            hasattr(value, '__getitem__') and
            hasattr(value, '__setitem__')
        )
        
        if is_maplike:
            from ..common import edit_maplike
            return edit_maplike(value, parent=parent)
        
        # Route generic array-like objects to NDarray editor
        has_array = hasattr(value, '__array__')
        if has_array:
            return NDarray.edit(value, parent=parent)
        
        # Fallback to Generic editor for other terminal types
        return Generic.edit(value, parent=parent)
    
    # For complex types with properties, use unified property editor
    ro_props, rw_props = inspect_type(value)
    
    if ro_props or rw_props:
        # Has properties - use unified editor
        return unified_edit(value, parent=parent)
    
    # No properties - try specialized editor or generic
    current_module = sys.modules[__name__]
    if hasattr(current_module, type_name):
        edit_module = getattr(current_module, type_name)
        if hasattr(edit_module, 'edit'):
            return edit_module.edit(value, parent=parent)
    
    # Final fallback to Generic
    return Generic.edit(value, parent=parent)
