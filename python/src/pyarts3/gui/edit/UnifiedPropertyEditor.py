"""Unified property-based editor for ARTS types.

This module provides a single-window, chainable editor that works by inspecting
object properties and presenting them in a scrollable list. Navigation behaves
like a file browser:

- A breadcrumb at the top shows the current path and lets you jump back
- The main area lists properties (read-write first, then read-only)
- Double-click a terminal value to edit it in-place (using existing editors)
- Double-click a complex value to navigate into it in the SAME window

This solves:
1) Too many properties to fit in tabs (now a scrollable table)
2) Navigation popping up new windows (now an in-dialog navigator with breadcrumbs)
"""

import numpy as np
from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QWidget,
    QPushButton, QDialogButtonBox, QLabel, QScrollArea,
    QTableWidget, QTableWidgetItem, QHeaderView, QApplication
)
from PyQt5.QtCore import Qt


# Properties to ignore when inspecting numpy arrays
_IGNORE_ARRAY = set(dir(np.ndarray))


def inspect_type(t):
    """Inspect an ARTS type and return read-only and read-write properties.
    
    Parameters
    ----------
    t : type
        The type to inspect
        
    Returns
    -------
    ro : list
        List of read-only property names
    rw : list
        List of read-write property names
    """
    import inspect as _inspect

    # Prefer class-level descriptor introspection to avoid triggering instance getters
    ro: list[str] = []
    rw: list[str] = []

    cls = type(t)
    for name in dir(cls):
        if name.startswith('_'):
            continue
        try:
            obj = getattr(cls, name)
        except Exception:
            continue

        # Skip callables (methods, functions)
        if _inspect.isroutine(obj) or isinstance(obj, (staticmethod, classmethod)):
            continue

        # Python @property
        if isinstance(obj, property):
            if getattr(obj, 'fset', None) is not None:
                rw.append(name)
            else:
                ro.append(name)
            continue

        # C-extension descriptors (e.g., nanobind/pybind11 getset/member descriptors)
        has_get = hasattr(obj, '__get__')
        has_set = hasattr(obj, '__set__')
        if has_get:
            (rw if has_set else ro).append(name)

    # If nothing found, do a conservative instance-level fallback:
    if not ro and not rw:
        try:
            t_members = dir(t)
        except Exception:
            t_members = []
        has_array = "__array__" in t_members

        for name in t_members:
            if name.startswith('_'):
                continue
            if has_array and name in _IGNORE_ARRAY:
                continue
            # Avoid calling heavy/callable attributes
            try:
                val = getattr(t, name)
            except Exception:
                continue
            # Skip callables
            try:
                if callable(val):
                    continue
            except Exception:
                pass
            # Treat as read-only by default (safer); editor will guard setting
            ro.append(name)

    return sorted(ro), sorted(rw)


def is_terminal_type(value):
    """Check if a value is a terminal type that can be edited directly.
    
    Terminal types are:
    - ARTS: Numeric, Index, String, Time, SpeciesIsotope
    - Python: int, float, bool, str
    - NumPy: integer, floating types
    - NumPy: ndarray
    - ARTS: ArrayOf types
    - ARTS: Map-like types (keys, items, values, __getitem__, __setitem__)
    - ARTS: Options enums
    - ARTS: GriddedField types
    
    Parameters
    ----------
    value : any
        The value to check
        
    Returns
    -------
    bool
        True if this is a terminal type
    """
    import pyarts3.arts as arts
    
    # Check type name for ARTS types
    type_name = type(value).__name__
    
    # ARTS terminal types
    if type_name in ('Numeric', 'Index', 'String', 'Time', 'SpeciesIsotope'):
        return True
    
    # ArrayOf types
    if type_name.startswith('ArrayOf'):
        return True
    
    # Check if it's a map-like type (has dict-like interface)
    try:
        is_maplike = (
            hasattr(value, 'keys') and
            hasattr(value, 'items') and
            hasattr(value, 'values') and
            hasattr(value, '__getitem__') and
            hasattr(value, '__setitem__')
        )
        if is_maplike:
            return True
    except Exception:
        pass
    
    # Check if it's a GriddedField type
    try:
        is_griddedfield = (
            hasattr(value, '__array__') and
            hasattr(value, 'grids') and
            hasattr(value, 'gridnames') and
            hasattr(value, 'dataname')
        )
        if is_griddedfield:
            return True
    except Exception:
        pass
    
    # Check if it's an option enum
    try:
        option_groups = arts.globals.option_groups()
        if type_name in option_groups:
            return True
    except Exception:
        pass
    
    # Python built-in types
    if isinstance(value, (int, float, bool, str)):
        return True
    
    # NumPy types
    if isinstance(value, (np.integer, np.floating)):
        return True
    
    # NumPy ndarray
    if isinstance(value, np.ndarray):
        return True
    
    # Check if it has __array__ (array-like) - but not if it's already
    # caught by GriddedField check above
    if hasattr(value, '__array__'):
        return True
    
    return False


def edit_terminal(value, parent=None, owner=None, prop_name=None):
    """Edit a terminal type using the appropriate specialized editor.
    
    This function routes terminal types to their specialized editors:
    - ARTS: Numeric, Index, String, Time, SpeciesIsotope → dedicated module editors
    - ArrayOf* → ArrayOf editor
    - Options enums → Options editor
    - GriddedField → edit_griddedfield
    - Map-like → edit_maplike
    - Array-like → edit_ndarraylike
    - Python primitives → Generic editor
    
    Parameters
    ----------
    value : any
        The terminal value to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    any or None
        Edited value if accepted, None if cancelled
    """
    import sys
    from . import Generic, ArrayOf, Options
    from . import Index, Numeric, String, Time, SpeciesIsotope
    
    type_name = type(value).__name__
    current_module = sys.modules['pyarts3.gui.edit']
    
    # Check if we have a specific editor module for this type
    if hasattr(current_module, type_name):
        edit_module = getattr(current_module, type_name)
        if hasattr(edit_module, 'edit'):
            return edit_module.edit(value, parent=parent)
    
    # Route ArrayOf* to the generic ArrayOf editor with type constraints
    if type_name.startswith('ArrayOf'):
        # Extract inner type from __doc__ (e.g., 'A list of :class:`~pyarts3.arts.Index`')
        allowed_item_types = None
        try:
            doc = type(value).__doc__
            if doc:
                # Look for :class:`~pyarts3.arts.TypeName` pattern
                import re
                match = re.search(r':class:`~pyarts3\.arts\.(\w+)`', doc)
                if match:
                    inner_type_name = match.group(1)
                    import pyarts3.arts as arts
                    if hasattr(arts, inner_type_name):
                        allowed_item_types = [getattr(arts, inner_type_name)]
        except Exception:
            pass
        
        # Use edit_listlike with type constraints instead of dedicated ArrayOf editor
        from .widgets import edit_listlike
        return edit_listlike(value, parent=parent, allowed_item_types=allowed_item_types)
    
    # Check if this is an option group enum
    try:
        import pyarts3.arts as arts
        option_groups = arts.globals.option_groups()
        if type_name in option_groups:
            return Options.edit(value, parent=parent)
    except Exception:
        pass
    
    # Check if this is a gridded field
    try:
        is_griddedfield = (
            hasattr(value, '__array__') and
            hasattr(value, 'grids') and
            hasattr(value, 'gridnames') and
            hasattr(value, 'dataname')
        )
        if is_griddedfield:
            from .widgets import edit_griddedfield
            return edit_griddedfield(value, parent=parent)
    except Exception:
        pass
    
    # Check if this is a map-like object
    try:
        is_maplike = (
            hasattr(value, 'keys') and
            hasattr(value, 'items') and
            hasattr(value, 'values') and
            hasattr(value, '__getitem__') and
            hasattr(value, '__setitem__')
        )
        if is_maplike:
            # Contract: Derive allowed types from:
            # 1. Property getter signature if owner/prop_name provided (dict[K, V] form)
            # 2. Direct type's __getitem__.__doc__ for builtin dict types
            # Otherwise read-only.
            def _resolve_type(name):
                try:
                    # Remove pyarts3.arts. prefix if present
                    if name.startswith('pyarts3.arts.'):
                        name = name[len('pyarts3.arts.'):]
                    # Try to resolve against pyarts3.arts namespace
                    import pyarts3.arts as arts
                    if hasattr(arts, name):
                        return getattr(arts, name)
                except Exception:
                    pass
                return None

            allowed_key_types = None
            allowed_value_types = None
            
            # Strategy 1: Parse property getter signature (for dict properties)
            try:
                if owner is not None and prop_name:
                    owner_cls = type(owner)
                    desc = getattr(owner_cls, prop_name, None)
                    if desc is not None:
                        fget = getattr(desc, 'fget', None)
                        if fget is not None and hasattr(fget, '__nb_signature__'):
                            sigs = getattr(fget, '__nb_signature__')
                            for entry in sigs:
                                sig = entry[0] if isinstance(entry, (list, tuple)) else str(entry)
                                # Look for '-> dict[KeyType, ValueType]'
                                arrow = sig.find('->')
                                if arrow == -1:
                                    continue
                                ret_type = sig[arrow + 2:].strip()
                                # Parse dict[K, V] form
                                if ret_type.startswith('dict['):
                                    inside = ret_type[5:]  # skip 'dict['
                                    if inside.endswith(']'):
                                        inside = inside[:-1]
                                    # Split by comma at top level
                                    parts = []
                                    buf = []
                                    depth = 0
                                    for ch in inside:
                                        if ch == '[':
                                            depth += 1
                                        elif ch == ']':
                                            depth -= 1
                                        if ch == ',' and depth == 0:
                                            parts.append(''.join(buf).strip())
                                            buf = []
                                        else:
                                            buf.append(ch)
                                    if buf:
                                        parts.append(''.join(buf).strip())
                                    if len(parts) >= 2:
                                        key_name = parts[0]
                                        val_name = parts[1]
                                        kt = _resolve_type(key_name)
                                        vt = _resolve_type(val_name)
                                        if kt:
                                            allowed_key_types = [kt]
                                        if vt:
                                            allowed_value_types = [vt]
                                        break
            except Exception:
                pass
            
            # Strategy 2: Parse __getitem__.__doc__ for builtin dict types
            if allowed_key_types is None or allowed_value_types is None:
                try:
                    t = type(value)
                    getitem = getattr(t, '__getitem__', None)
                    if getitem is not None:
                        doc = getattr(getitem, '__doc__', None)
                        if doc:
                            # Parse: __getitem__(self, arg: KeyType, /) -> ValueType
                            arrow = doc.find('->')
                            if arrow != -1:
                                # Extract value type from return
                                ret_part = doc[arrow + 2:].strip()
                                # Remove trailing period if present
                                if ret_part.endswith('.'):
                                    ret_part = ret_part[:-1]
                                vt = _resolve_type(ret_part)
                                if vt and allowed_value_types is None:
                                    allowed_value_types = [vt]
                            
                            # Extract key type from arg:
                            arg_start = doc.find('arg:')
                            if arg_start != -1:
                                arg_part = doc[arg_start + 4:].strip()
                                # Find end (comma or slash)
                                for end_char in [',', '/']:
                                    end_idx = arg_part.find(end_char)
                                    if end_idx != -1:
                                        arg_part = arg_part[:end_idx].strip()
                                        break
                                kt = _resolve_type(arg_part)
                                if kt and allowed_key_types is None:
                                    allowed_key_types = [kt]
                except Exception:
                    pass

            from .widgets import edit_maplike
            return edit_maplike(
                value,
                parent=parent,
                allowed_key_types=allowed_key_types,
                allowed_value_types=allowed_value_types,
            )
    except Exception:
        pass
    
    # Route generic array-like objects to ndarraylike editor
    try:
        if hasattr(value, '__array__'):
            from .widgets import edit_ndarraylike
            # Attempt to determine allowed_shape for ARTS-derived arrays.
            # Strategy: default-construct the type and read its shape; per design,
            # zeros indicate dynamic dims and non-zeros are fixed.
            allowed_shape = None
            try:
                t = type(value)
                # Avoid numpy.ndarray default construction which is not supported
                is_numpy_arr = (t is getattr(np, 'ndarray', None))
            except Exception:
                is_numpy_arr = False
            if not is_numpy_arr:
                try:
                    default_inst = type(value)()
                    shp = getattr(default_inst, 'shape', None)
                    if shp is None:
                        try:
                            shp = np.array(default_inst).shape
                        except Exception:
                            shp = None
                    if shp and isinstance(shp, tuple) and len(shp) == np.array(value).ndim:
                        # Normalize to ints
                        allowed_shape = tuple(int(x) for x in shp)
                except Exception:
                    allowed_shape = None
            return edit_ndarraylike(value, parent=parent, allowed_shape=allowed_shape)
    except Exception:
        pass
    
    # Fallback to Generic editor for other terminal types
    return Generic.edit(value, parent=parent)


class PropertyEditor(QDialog):
    """Unified property-based editor for ARTS types.
    
    Single-window navigator:
    - Breadcrumb at top with clickable segments
    - Scrollable list of properties (RW first, then RO)
    - Double-click to edit terminal or navigate into complex types
    """

    def __init__(self, value, parent=None):
        super().__init__(parent)
        self.setMinimumSize(900, 650)

        # Root working copy (we never mutate the original directly)
        self._root = self._copy_value(value)
        self._type_name = type(value).__name__

        # Navigation state
        # path: list[str] of property names from root to current
        self._path: list[str] = []
        # readonly flags along the path (propagates RO to children)
        self._ro_flags: list[bool] = []

        self.setWindowTitle(f"Edit {self._type_name}")

        # UI
        self._outer = QVBoxLayout(self)

        # Breadcrumb bar
        self._crumb_bar = QHBoxLayout()
        self._outer.addLayout(self._crumb_bar)

        # Property table
        self._table = QTableWidget(0, 4)
        self._table.setHorizontalHeaderLabels(["Name", "Type", "Access", "Value preview"])
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self._table.setSelectionBehavior(self._table.SelectRows)
        self._table.setEditTriggers(QTableWidget.NoEditTriggers)
        self._table.cellDoubleClicked.connect(self._on_double_click)
        self._outer.addWidget(self._table)

        # Dialog buttons
        self._buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self._buttons.accepted.connect(self.accept)
        self._buttons.rejected.connect(self.reject)
        self._outer.addWidget(self._buttons)

        # Initial render
        self._refresh_view()

    def _copy_value(self, value):
        import copy
        try:
            return copy.deepcopy(value)
        except Exception:
            try:
                return copy.copy(value)
            except Exception:
                return value

    # -------- Navigation helpers --------
    def _resolve_current(self):
        """Resolve the current object by walking _path from _root."""
        obj = self._root
        for name in self._path:
            try:
                obj = getattr(obj, name)
            except Exception:
                break
        return obj

    def _current_readonly(self) -> bool:
        return any(self._ro_flags) if self._ro_flags else False

    def _push_path(self, prop_name: str, parent_rw_props: list[str]):
        is_ro = prop_name not in parent_rw_props
        self._path.append(prop_name)
        self._ro_flags.append(is_ro or self._current_readonly())
        self._refresh_view()

    def _pop_to(self, depth: int):
        """Pop navigation to given depth (0 = root)."""
        if depth < 0:
            depth = 0
        self._path = self._path[:depth]
        self._ro_flags = self._ro_flags[:depth]
        self._refresh_view()

    # -------- UI refresh --------
    def _refresh_view(self):
        # Breadcrumbs
        # Clear layout
        while self._crumb_bar.count():
            item = self._crumb_bar.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

        # Root button
        root_btn = QPushButton(self._type_name)
        root_btn.setFlat(True)
        root_btn.clicked.connect(lambda: self._pop_to(0))
        self._crumb_bar.addWidget(root_btn)

        # Add separators and buttons for each path segment
        for idx, name in enumerate(self._path):
            sep = QLabel(" > ")
            self._crumb_bar.addWidget(sep)
            btn = QPushButton(name)
            btn.setFlat(True)
            btn.clicked.connect(lambda _, d=idx + 1: self._pop_to(d))
            self._crumb_bar.addWidget(btn)

        self._crumb_bar.addStretch()

        # Table of properties
        current_obj = self._resolve_current()
        ro_props, rw_props = inspect_type(current_obj)

        # Build ordered list: RW first, then RO
        ordered = [(name, False) for name in rw_props] + [(name, True) for name in ro_props]

        self._table.setRowCount(len(ordered))
        for row, (name, is_ro) in enumerate(ordered):
            # Name
            name_item = QTableWidgetItem(name + (" 🔒" if is_ro else ""))
            self._table.setItem(row, 0, name_item)

            # Type and value
            try:
                val = getattr(current_obj, name)
                tname = type(val).__name__
                preview = self._preview_value(val)
            except Exception as e:
                val = None
                tname = "<error>"
                preview = str(e)

            type_item = QTableWidgetItem(tname)
            self._table.setItem(row, 1, type_item)

            # Access
            access_item = QTableWidgetItem("RO" if (is_ro or self._current_readonly()) else "RW")
            self._table.setItem(row, 2, access_item)

            # Preview
            preview_item = QTableWidgetItem(preview)
            self._table.setItem(row, 3, preview_item)

        self._table.resizeRowsToContents()

    def _preview_value(self, v, maxlen: int = 120) -> str:
        try:
            if hasattr(v, '__array__') and not isinstance(v, (str, bytes)):
                # Short preview for arrays
                try:
                    shape = getattr(v, 'shape', None)
                    if shape is None and hasattr(v, '__array__'):
                        arr = np.array(v)
                        shape = arr.shape
                    return f"array-like shape={shape}"
                except Exception:
                    return f"array-like"
            s = str(v)
        except Exception as e:
            s = f"<error: {e}>"
        if len(s) > maxlen:
            s = s[:maxlen - 3] + "..."
        return s

    # -------- Editing / Navigation --------
    def _on_double_click(self, row: int, col: int):
        current_obj = self._resolve_current()
        ro_props, rw_props = inspect_type(current_obj)
        ordered = [(name, False) for name in rw_props] + [(name, True) for name in ro_props]
        if row < 0 or row >= len(ordered):
            return
        name, is_ro = ordered[row]

        try:
            val = getattr(current_obj, name)
        except Exception:
            return

        # Determine access considering inherited RO status
        effective_ro = is_ro or self._current_readonly()

        # Terminal types → edit directly (if not RO), else view-only
        if is_terminal_type(val):
            if effective_ro:
                # View-only: still allow opening the editor but do not set back
                _ = edit_terminal(val, parent=self, owner=current_obj, prop_name=name)
                return
            new_val = edit_terminal(val, parent=self, owner=current_obj, prop_name=name)
            if new_val is not None:
                self._set_current_property(name, new_val)
                self._refresh_view()
            return

        # Complex types → navigate into same window
        # Push into path; RO propagates from parent RW list
        self._push_path(name, rw_props)

    def _set_current_property(self, prop_name: str, new_value):
        """Set a property on the current object in the navigation path.

        Always set on the CURRENT object, not the parent. This avoids collisions
        when parent and child share property names (e.g., AbsorptionLine.gl vs
        ZeemanLineModel.gl when editing 'z.gl').
        """
        target_obj = self._resolve_current()
        try:
            setattr(target_obj, prop_name, new_value)
        except Exception:
            # Silently ignore errors; some properties may be read-only or have validation
            pass

    def get_result(self):
        return self._root


def edit(value, parent=None):
    """Unified edit function for ARTS types.
    
    This function provides a chainable, property-based editing interface.
    It automatically determines whether to use terminal editing or create
    a navigable property explorer.
    
    Parameters
    ----------
    value : any
        The value to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    any or None
        Edited value if accepted, None if cancelled
    """
    # Ensure a QApplication exists to avoid Qt aborts in scripts/REPL
    app = QApplication.instance()
    if app is None:
        app = QApplication([])

    # Guard: Workspace must not be handled by unified; delegate to dedicated editor
    if type(value).__name__ == 'Workspace':
        try:
            from . import Workspace as WorkspaceEditor
            return WorkspaceEditor.edit(value, parent=parent)
        except Exception as e:
            # Fallback: raise a clear error instead of risking a segfault
            raise RuntimeError("Workspace cannot be edited with UnifiedPropertyEditor; use Workspace editor") from e

    # Check if terminal type
    if is_terminal_type(value):
        return edit_terminal(value, parent=parent)
    
    # Check if it has properties to edit
    ro_props, rw_props = inspect_type(value)

    # Always open the property editor for non-terminal types, even if we
    # couldn't detect properties (some types hide descriptors); the dialog
    # will render an empty list safely instead of misrouting to terminal.
    editor = PropertyEditor(value, parent=parent)
    if editor.exec_() == QDialog.Accepted:
        return editor.get_result()
    
    return None
