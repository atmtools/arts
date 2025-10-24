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
    QTableWidget, QTableWidgetItem, QHeaderView
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
    - ARTS: Numeric, Index, String
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
    if type_name in ('Numeric', 'Index', 'String'):
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


def edit_terminal(value, parent=None):
    """Edit a terminal type using the existing specialized editors.
    
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
    # Import the existing dispatch system
    from . import edit as dispatch_edit
    return dispatch_edit(value, parent=parent)


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
                _ = edit_terminal(val, parent=self)
                return
            new_val = edit_terminal(val, parent=self)
            if new_val is not None:
                self._set_current_property(name, new_val)
                self._refresh_view()
            return

        # Complex types → navigate into same window
        # Push into path; RO propagates from parent RW list
        self._push_path(name, rw_props)

    def _set_current_property(self, prop_name: str, new_value):
        """Set property on the object at current path (parent level)."""
        # Resolve parent object (one level up)
        parent_obj = self._root
        if len(self._path) > 0:
            for name in self._path[:-1]:
                parent_obj = getattr(parent_obj, name)
        target_obj = self._resolve_current() if len(self._path) == 0 else parent_obj

        try:
            setattr(target_obj, prop_name, new_value)
        except Exception as e:
            # Fallback: if setting failed on parent/target logic, try direct
            try:
                setattr(self._resolve_current(), prop_name, new_value)
            except Exception:
                print(f"Warning: Could not set {prop_name}: {e}")

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
