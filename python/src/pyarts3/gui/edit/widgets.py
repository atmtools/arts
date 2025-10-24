"""Editor widgets and helper functions for ARTS types.

This module contains all the editor functions and widgets that were previously
scattered in common.py. These are specific to editing ARTS types.
"""

import numpy as np
from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QPushButton,
    QDialogButtonBox, QSpinBox, QDoubleSpinBox, QLineEdit, QTextEdit, QWidget,
    QTableWidget, QTableWidgetItem, QHeaderView, QComboBox, QToolBar, QAction, QMessageBox
)
from PyQt5.QtCore import Qt


__all__ = [
    'ScientificDoubleSpinBox',
    'create_numeric_spinbox',
    'edit_numeric',
    'edit_index',
    'edit_string',
    'edit_listlike',
    'edit_ndarraylike',
    'edit_griddedfield',
    'edit_maplike',
    'create_atm_keys_table',
    'create_subsurface_keys_table',
    'create_surface_keys_table',
]


class ScientificDoubleSpinBox(QDoubleSpinBox):
    """QDoubleSpinBox that accepts scientific notation in text input.
    
    This widget provides the best of both worlds:
    - Arrow buttons for incrementing/decrementing
    - Text input that accepts intermediate invalid states like "1e" while typing "1e3"
    
    Unlike the default QDoubleSpinBox, this allows typing scientific notation
    naturally without validation errors.
    """
    
    def __init__(self, parent=None):
        super().__init__(parent)
        from PyQt5.QtGui import QDoubleValidator
        from PyQt5.QtCore import QLocale
        
        # Configure for large range
        self.setRange(-1e308, 1e308)
        self.setDecimals(10)
        self.setKeyboardTracking(False)
        
        # Replace the default validator with one that accepts scientific notation
        # and intermediate states
        validator = QDoubleValidator()
        validator.setNotation(QDoubleValidator.ScientificNotation)
        validator.setLocale(QLocale.c())
        self.lineEdit().setValidator(validator)
    
    def textFromValue(self, value):
        """Format value for display - use scientific notation for large/small values."""
        if abs(value) >= 1e6 or (abs(value) < 1e-3 and value != 0):
            return f"{value:.6e}"
        else:
            return f"{value:.10g}"
    
    def valueFromText(self, text):
        """Parse text to value - supports scientific notation."""
        text = text.strip()
        if not text:
            return 0.0
        try:
            return float(text)
        except ValueError:
            # If parsing fails, return current value
            return self.value()


def create_numeric_spinbox(value=0.0):
    """
    Create a configured spinbox for numeric editing with scientific notation support.
    
    Parameters
    ----------
    value : float
        Initial value
    
    Returns
    -------
    ScientificDoubleSpinBox
        Configured spinbox widget with arrows and flexible text input
    """
    spin_box = ScientificDoubleSpinBox()
    spin_box.setValue(float(value))
    return spin_box


def edit_numeric(value, parent=None):
    """
    Edit a numeric value with arrows and flexible text input.
    
    Provides a spinbox that:
    - Has arrow buttons for incrementing/decrementing
    - Accepts scientific notation in text field (e.g., 1e3, 2.5e-10)
    - Allows intermediate invalid states while typing (e.g., "1e" while typing "1e3")
    
    Parameters
    ----------
    value : float or Numeric
        The value to edit
    parent : QWidget, optional
        Parent widget
    
    Returns
    -------
    float or None
        Edited value if accepted, None if cancelled
    """
    from .common import create_simple_editor_dialog
    
    spin_box = ScientificDoubleSpinBox()
    spin_box.setValue(float(value))
    
    return create_simple_editor_dialog(
        "Edit Numeric",
        spin_box,
        spin_box.value,
        parent
    )


def edit_index(value, parent=None):
    """
    Edit an integer value using a spin box.
    
    Parameters
    ----------
    value : int or Index
        The value to edit
    parent : QWidget, optional
        Parent widget
    
    Returns
    -------
    int or None
        Edited value if accepted, None if cancelled
    """
    from .common import create_simple_editor_dialog
    
    spin_box = QSpinBox()
    spin_box.setRange(-2147483648, 2147483647)
    spin_box.setValue(int(value))
    
    return create_simple_editor_dialog(
        "Edit Index",
        spin_box,
        spin_box.value,
        parent
    )


def edit_string(value, parent=None):
    """
    Edit a string value using a line edit.
    
    Parameters
    ----------
    value : str or String
        The value to edit
    parent : QWidget, optional
        Parent widget
    
    Returns
    -------
    str or None
        Edited value if accepted, None if cancelled
    """
    from .common import create_simple_editor_dialog
    
    line_edit = QLineEdit()
    line_edit.setText(str(value))
    
    return create_simple_editor_dialog(
        "Edit String",
        line_edit,
        line_edit.text,
        parent
    )


def edit_listlike(value, parent=None):
    """
    Edit a list-like value (list, tuple, etc.) in a simple table.
    Shows index and value columns with inline editing.
    
    Parameters
    ----------
    value : list or tuple
        The list-like value to edit
    parent : QWidget, optional
        Parent widget
    
    Returns
    -------
    list or None
        Edited value as list if accepted, None if cancelled
    """
    from PyQt5.QtWidgets import QApplication
    
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit List")
    
    layout = QVBoxLayout()
    
    # Convert to list for editing
    data = list(value)
    
    # Table
    table = QTableWidget()
    table.setColumnCount(2)
    table.setHorizontalHeaderLabels(["Index", "Value"])
    table.horizontalHeader().setStretchLastSection(True)
    
    table.setRowCount(len(data))
    for i, val in enumerate(data):
        # Index (read-only)
        idx_item = QTableWidgetItem(str(i))
        idx_item.setFlags(idx_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(i, 0, idx_item)
        
        # Value (editable) - try to format nicely
        if isinstance(val, (int, float)):
            val_item = QTableWidgetItem(f"{val:.10g}" if isinstance(val, float) else str(val))
        else:
            val_item = QTableWidgetItem(str(val))
        table.setItem(i, 1, val_item)
    
    # Add double-click handler for numeric editing if all values are numeric
    all_numeric = all(isinstance(v, (int, float, np.number)) for v in data)
    if all_numeric:
        def on_cell_double_clicked(row, col):
            if col == 1:  # Only value column
                try:
                    current_val = float(table.item(row, 1).text())
                    new_val = edit_numeric(current_val, dialog)
                    if new_val is not None:
                        table.item(row, 1).setText(f"{new_val:.10g}")
                except (ValueError, TypeError):
                    pass  # Skip if not numeric
        
        table.cellDoubleClicked.connect(on_cell_double_clicked)
    
    layout.addWidget(table)
    
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    if dialog.exec_() == QDialog.Accepted:
        result = []
        for i in range(table.rowCount()):
            text = table.item(i, 1).text()
            # Try to preserve type
            if all_numeric:
                try:
                    # Try int first, then float
                    if '.' not in text and 'e' not in text.lower():
                        result.append(int(text))
                    else:
                        result.append(float(text))
                except (ValueError, AttributeError):
                    result.append(data[i])  # Keep original on error
            else:
                result.append(text)
        return result
    return None


def edit_ndarraylike(value, parent=None):
    """
    Edit an ndarray-like value (numpy array, ARTS arrays, etc.) in a table.
    Handles 1D and 2D arrays.
    
    Parameters
    ----------
    value : ndarray or array-like
        The array value to edit
    parent : QWidget, optional
        Parent widget
    
    Returns
    -------
    numpy.ndarray or None
        Edited value as numpy array if accepted, None if cancelled
    """
    from PyQt5.QtWidgets import QApplication
    
    # Ensure a QApplication exists (avoid Qt crashes if none is running)
    app = QApplication.instance()
    if app is None:
        # Create an application with no argv to minimize side effects
        app = QApplication([])
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit Array")
    
    layout = QVBoxLayout()
    
    # Capture original type for reconstruction
    original_type = type(value)
    
    # Convert to numpy array if needed
    if hasattr(value, '__array__'):
        data = np.array(value)
    else:
        data = np.array(value)
    
    # Handle 0D arrays (scalars) - dispatch back to edit() to preserve type
    if data.ndim == 0:
        # Import here to avoid circular dependency
        from . import edit as edit_module
        # Extract scalar value (preserves int/float type)
        scalar_value = data.item()
        result = edit_module.edit(scalar_value, parent=parent)
        if result is None:
            return None
        # Reconstruct original type if it changed
        original_type = type(value)
        result_type = type(result)
        if original_type != result_type:
            try:
                return original_type(result)
            except Exception:
                # If conversion fails, return as-is
                return result
        return result
    
    # Handle 1D arrays
    if data.ndim == 1:
        dialog.resize(600, 400)
        
        # Info label
        info = QLabel(f"Size: {len(data)}")
        layout.addWidget(info)
        
        # Table
        table = QTableWidget()
        table.setColumnCount(2)
        table.setHorizontalHeaderLabels(["Index", "Value"])
        table.horizontalHeader().setStretchLastSection(True)
        
        table.setRowCount(len(data))
        for i, val in enumerate(data):
            # Index (read-only)
            idx_item = QTableWidgetItem(str(i))
            idx_item.setFlags(idx_item.flags() & ~Qt.ItemIsEditable)
            table.setItem(i, 0, idx_item)
            
            # Value (editable) - store original value in UserRole
            val_item = QTableWidgetItem(f"{val:.10g}")
            val_item.setData(Qt.UserRole, float(val))  # Store original
            table.setItem(i, 1, val_item)
        
        # Track which cells have been edited
        edited_cells = set()
        
        # Add double-click handler for numeric editing
        def on_cell_double_clicked(row, col):
            if col == 1:  # Only value column
                current_val = table.item(row, 1).data(Qt.UserRole)
                new_val = edit_numeric(current_val, dialog)
                if new_val is not None:
                    table.item(row, 1).setText(f"{new_val:.10g}")
                    table.item(row, 1).setData(Qt.UserRole, float(new_val))
                    edited_cells.add(row)
        
        table.cellDoubleClicked.connect(on_cell_double_clicked)
        
        layout.addWidget(table)
        
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        
        dialog.setLayout(layout)
        
        if dialog.exec_() == QDialog.Accepted:
            result = []
            for i in range(table.rowCount()):
                item = table.item(i, 1)
                # Use stored original value (which gets updated on edit)
                result.append(item.data(Qt.UserRole))
            result_array = np.array(result)
            # Reconstruct original type if needed
            if original_type != np.ndarray:
                try:
                    return original_type(result_array)
                except Exception:
                    return result_array
            return result_array
        return None
    
    # Handle 2D arrays
    elif data.ndim == 2:
        dialog.resize(800, 600)
        
        # Handle empty arrays
        if data.size == 0:
            layout.addWidget(QLabel(f"Empty array with shape: {data.shape}"))
            layout.addWidget(QLabel("Cannot edit empty arrays."))
            layout.addWidget(QLabel("Please reshape or populate the array first."))
            buttons = QDialogButtonBox(QDialogButtonBox.Ok)
            buttons.accepted.connect(dialog.accept)
            layout.addWidget(buttons)
            dialog.setLayout(layout)
            dialog.exec_()
            return None
        
        # Info label
        info = QLabel(f"Shape: {data.shape[0]} × {data.shape[1]}")
        layout.addWidget(info)
        
        # Table
        table = QTableWidget()
        table.setRowCount(data.shape[0])
        table.setColumnCount(data.shape[1])
        
        # Set headers
        table.setHorizontalHeaderLabels([str(i) for i in range(data.shape[1])])
        table.setVerticalHeaderLabels([str(i) for i in range(data.shape[0])])
        
        # Fill data - store original values in UserRole
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                item = QTableWidgetItem(f"{data[i, j]:.6g}")
                item.setData(Qt.UserRole, float(data[i, j]))  # Store original
                table.setItem(i, j, item)
        
        # Add double-click handler for numeric editing
        def on_cell_double_clicked(row, col):
            item = table.item(row, col)
            current_val = item.data(Qt.UserRole)
            new_val = edit_numeric(current_val, dialog)
            if new_val is not None:
                item.setText(f"{new_val:.6g}")
                item.setData(Qt.UserRole, float(new_val))
        
        table.cellDoubleClicked.connect(on_cell_double_clicked)
        
        layout.addWidget(table)
        
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        
        dialog.setLayout(layout)
        
        if dialog.exec_() == QDialog.Accepted:
            # Build result array from stored values
            result = []
            for i in range(data.shape[0]):
                row = []
                for j in range(data.shape[1]):
                    item = table.item(i, j)
                    row.append(item.data(Qt.UserRole))
                result.append(row)
            result_array = np.array(result)
            # Reconstruct original type if needed
            if original_type != np.ndarray:
                try:
                    return original_type(result_array)
                except Exception:
                    return result_array
            return result_array
        return None
    
    # Handle higher dimensions (3D+) - show error
    else:
        layout.addWidget(QLabel(f"Array has {data.ndim} dimensions"))
        layout.addWidget(QLabel("Can only edit 1D and 2D arrays"))
        layout.addWidget(QLabel("Please reshape the array first"))
        buttons = QDialogButtonBox(QDialogButtonBox.Ok)
        buttons.accepted.connect(dialog.accept)
        layout.addWidget(buttons)
        dialog.setLayout(layout)
        dialog.exec_()
        return None


def edit_griddedfield(value, parent=None):
    """
    Edit a GriddedField value (grids, gridnames, dataname, data).
    
    Parameters
    ----------
    value : GriddedField*
        Any GriddedField type
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    same-type-as-input or None
        The edited GriddedField if accepted, None if cancelled
    """
    from PyQt5.QtWidgets import QApplication
    
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    
    original_type = type(value)
    
    dialog = QDialog(parent)
    dialog.setWindowTitle(f"Edit {original_type.__name__}")
    dialog.resize(900, 700)
    
    layout = QVBoxLayout()
    
    # Info
    info = QLabel(f"Gridded field with {len(value.grids)} grid(s)")
    info.setStyleSheet("color: #555; font-weight: bold;")
    layout.addWidget(info)
    
    # Check if grids/gridnames are tuples (some GriddedField types use tuples, some lists)
    grids_was_tuple = isinstance(value.grids, tuple)
    gridnames_was_tuple = isinstance(value.gridnames, tuple)
    
    # Convert to lists for editing
    edited_grids = list(value.grids) if value.grids else []
    edited_gridnames = list(value.gridnames) if value.gridnames else []
    
    # Grids table
    grids_label = QLabel("Grids (double-click to edit):")
    layout.addWidget(grids_label)
    
    grids_table = QTableWidget()
    grids_table.setColumnCount(3)
    grids_table.setHorizontalHeaderLabels(["Index", "Name", "Grid"])
    grids_table.horizontalHeader().setStretchLastSection(True)
    
    def refresh_grids_table():
        grids_table.setRowCount(len(edited_grids))
        for i, (grid, name) in enumerate(zip(edited_grids, edited_gridnames)):
            # Index
            idx_item = QTableWidgetItem(str(i))
            idx_item.setFlags(idx_item.flags() & ~Qt.ItemIsEditable)
            grids_table.setItem(i, 0, idx_item)
            
            # Name
            name_item = QTableWidgetItem(str(name))
            name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable)
            grids_table.setItem(i, 1, name_item)
            
            # Grid preview
            try:
                grid_arr = np.array(grid)
                grid_preview = f"shape={grid_arr.shape}, min={grid_arr.min():.3g}, max={grid_arr.max():.3g}"
            except Exception:
                grid_preview = f"<{type(grid).__name__}>"
            grid_item = QTableWidgetItem(grid_preview)
            grid_item.setData(Qt.UserRole, grid)
            grid_item.setFlags(grid_item.flags() & ~Qt.ItemIsEditable)
            grids_table.setItem(i, 2, grid_item)
    
    def on_grid_double_clicked(row, col):
        if row >= len(edited_grids):
            return
        current_grid = edited_grids[row]
        # Use edit_ndarraylike directly to avoid circular import issues
        new_grid = edit_ndarraylike(current_grid, parent=dialog)
        if new_grid is not None:
            edited_grids[row] = new_grid
            refresh_grids_table()
    
    grids_table.cellDoubleClicked.connect(on_grid_double_clicked)
    refresh_grids_table()
    layout.addWidget(grids_table)
    
    # Dataname field
    dataname_label = QLabel("Data name:")
    layout.addWidget(dataname_label)
    dataname_edit = QLineEdit()
    dataname_edit.setText(str(value.dataname) if value.dataname else "")
    layout.addWidget(dataname_edit)
    
    # Data table (read-only summary with edit button)
    data_label = QLabel("Data (double-click to edit):")
    layout.addWidget(data_label)
    
    edited_data = [value.data]  # Store in list so we can mutate it
    
    data_table = QTableWidget()
    data_table.setColumnCount(3)
    data_table.setHorizontalHeaderLabels(["Type", "Shape", "Summary"])
    data_table.horizontalHeader().setStretchLastSection(True)
    data_table.setRowCount(1)
    
    def refresh_data_table():
        data = edited_data[0]
        # Type
        type_item = QTableWidgetItem(type(data).__name__)
        type_item.setFlags(type_item.flags() & ~Qt.ItemIsEditable)
        data_table.setItem(0, 0, type_item)
        
        # Shape
        shape_item = QTableWidgetItem(str(data.shape))
        shape_item.setFlags(shape_item.flags() & ~Qt.ItemIsEditable)
        data_table.setItem(0, 1, shape_item)
        
        # Summary
        data_array = np.array(data)
        try:
            summary = f"min={data_array.min():.3g}, max={data_array.max():.3g}, mean={data_array.mean():.3g}"
        except Exception:
            summary = "(uninitialized)"
        summary_item = QTableWidgetItem(summary)
        summary_item.setData(Qt.UserRole, data)
        summary_item.setFlags(summary_item.flags() & ~Qt.ItemIsEditable)
        data_table.setItem(0, 2, summary_item)
    
    def on_data_double_clicked(row, col):
        item = data_table.item(0, 2)
        if item is None:
            return
        current_data = item.data(Qt.UserRole)
        # Use edit_ndarraylike directly
        new_data = edit_ndarraylike(current_data, parent=dialog)
        if new_data is not None:
            edited_data[0] = new_data
            refresh_data_table()
    
    data_table.cellDoubleClicked.connect(on_data_double_clicked)
    refresh_data_table()
    layout.addWidget(data_table)
    
    # Dialog buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    if dialog.exec_() == QDialog.Accepted:
        # Build result gridded field
        try:
            result = original_type()
            # Convert back to original container type if needed
            result.grids = tuple(edited_grids) if grids_was_tuple else edited_grids
            result.gridnames = tuple(edited_gridnames) if gridnames_was_tuple else edited_gridnames
            result.dataname = dataname_edit.text()
            result.data = edited_data[0]
            return result
        except Exception:
            return None
    return None


def edit_maplike(value, parent=None):
    """
    Edit a map-like value (dict-like with keys(), items(), values()).
    
    Map-like types have:
    - keys(): iterator over keys
    - items(): iterator over (key, value) pairs
    - values(): iterator over values
    - __getitem__, __setitem__, __delitem__
    
    Parameters
    ----------
    value : map-like instance
        Any ARTS map type (e.g., QuantumIdentifierNumericMap, WsvMap, etc.)
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    same-type-as-input or None
        The edited map if accepted, None if cancelled
    """
    from PyQt5.QtWidgets import QApplication
    from . import edit as dispatch_edit
    
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    
    def _parse_nb_signature_params(sig_obj):
        """
        Try to extract parameter type names from a nanobind __nb_signature__ entry.
        Returns a list of strings for parameter types (excluding return type).
        """
        try:
            # Expected shape from tests: iterable of tuples, first item is signature string
            sig = sig_obj[0] if isinstance(sig_obj, (list, tuple)) else str(sig_obj)
            # Find parameter list between the first '(' and the last ')'
            l = sig.find('(')
            r = sig.rfind(')')
            if l == -1 or r == -1 or r <= l:
                return []
            params_str = sig[l+1:r].strip()
            if not params_str:
                return []
            # Split parameters by comma, strip qualifiers like const&, &&, & and default names
            params = []
            for p in params_str.split(','):
                p = p.strip()
                # Drop common qualifiers
                p = p.replace('const', '').replace('&', '').replace('&&', '').strip()
                # Remove possible parameter names after a space
                # e.g., "Index i" -> "Index"
                if ' ' in p:
                    p = p.split(' ')[0]
                params.append(p)
            return [p for p in params if p]
        except Exception:
            return []
    
    def _resolve_type_name(name):
        """Resolve a C++/nanobind-exposed type name to a Python class or builtin."""
        try:
            # Common builtins mapping
            if name in ('int', 'size_t', 'Index'):
                # ARTS Index is exposed, try that first then fallback to int
                try:
                    import pyarts3.arts as arts
                    return getattr(arts, 'Index')
                except Exception:
                    return int
            if name in ('double', 'float', 'Numeric'):
                try:
                    import pyarts3.arts as arts
                    return getattr(arts, 'Numeric')
                except Exception:
                    return float
            if name in ('std::string', 'String', 'string', 'str'):
                try:
                    import pyarts3.arts as arts
                    return getattr(arts, 'String')
                except Exception:
                    return str
            # Try to resolve against pyarts3.arts namespace
            import pyarts3.arts as arts
            if hasattr(arts, name):
                return getattr(arts, name)
        except Exception:
            pass
        return None
    
    def _infer_map_types(m):
        """
        Infer possible (key_types, value_types) from __getitem__/__setitem__ nanobind signatures.
        Returns two lists of Python classes (may be empty).
        """
        key_types = []
        val_types = []
        try:
            getitem = getattr(type(m), '__getitem__', None) or getattr(m, '__getitem__', None)
            if getitem is not None and hasattr(getitem, '__nb_signature__'):
                for entry in getattr(getitem, '__nb_signature__'):
                    params = _parse_nb_signature_params(entry)
                    # Ignore 'self' if present and extract first parameter as key type
                    if len(params) >= 1:
                        t = _resolve_type_name(params[0])
                        if t and t not in key_types:
                            key_types.append(t)
        except Exception:
            pass
        try:
            setitem = getattr(type(m), '__setitem__', None) or getattr(m, '__setitem__', None)
            if setitem is not None and hasattr(setitem, '__nb_signature__'):
                for entry in getattr(setitem, '__nb_signature__'):
                    params = _parse_nb_signature_params(entry)
                    # Expect (key, value) parameters
                    if len(params) >= 2:
                        t = _resolve_type_name(params[1])
                        if t and t not in val_types:
                            val_types.append(t)
        except Exception:
            pass
        return key_types, val_types
    
    def _choose_type_dialog(type_list, title, parent_dialog):
        """Prompt the user to choose a type from a list. Returns the chosen class or None."""
        if not type_list:
            return None
        if len(type_list) == 1:
            return type_list[0]
        # Build a small modal dialog with a combo box
        dlg = QDialog(parent_dialog)
        dlg.setWindowTitle(title)
        v = QVBoxLayout(dlg)
        v.addWidget(QLabel("Select a type:"))
        combo = QComboBox(dlg)
        names = [t.__name__ if hasattr(t, '__name__') else str(t) for t in type_list]
        combo.addItems(names)
        v.addWidget(combo)
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        v.addWidget(buttons)
        if dlg.exec_() == QDialog.Accepted:
            return type_list[combo.currentIndex()]
        return None
    
    def _construct_default(t):
        """Try to construct a sensible default instance for type t."""
        try:
            return t()
        except Exception:
            # Fallback for some simple exposed types
            if t in (int, float, str):
                return t()
            # Try Numeric/Index explicit constructors
            try:
                import pyarts3.arts as arts
                if t is getattr(arts, 'Numeric', None):
                    return t(0.0)
                if t is getattr(arts, 'Index', None):
                    return t(0)
            except Exception:
                pass
        return None
    
    original_type = type(value)
    
    dialog = QDialog(parent)
    dialog.setWindowTitle(f"Edit {original_type.__name__}")
    dialog.resize(900, 600)
    
    layout = QVBoxLayout()
    
    # Info label
    info = QLabel(f"Map with {len(value)} entries")
    info.setStyleSheet("color: #555; font-weight: bold;")
    layout.addWidget(info)
    
    # Store edited entries
    # Convert to list of (key, value) pairs for editing
    edited_entries = list(value.items()) if len(value) > 0 else []
    deleted_indices = set()
    
    # Table for entries
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Key", "Value", "Types"])
    table.horizontalHeader().setStretchLastSection(True)
    table.setSelectionBehavior(QTableWidget.SelectRows)
    table.setSelectionMode(QTableWidget.SingleSelection)
    
    def _summarize(obj, max_len=80):
        try:
            s = str(obj)
        except Exception:
            s = f"<{type(obj).__name__}>"
        if len(s) > max_len:
            s = s[:max_len-1] + "…"
        return s
    
    def refresh_table():
        # Filter out deleted entries
        active_entries = [e for i, e in enumerate(edited_entries) if i not in deleted_indices]
        
        n = len(active_entries)
        info.setText(f"Map with {n} entries")
        
        table.setRowCount(n)
        for i, (key, val) in enumerate(active_entries):
            # Key
            key_item = QTableWidgetItem(_summarize(key))
            key_item.setData(Qt.UserRole, key)
            key_item.setFlags(key_item.flags() & ~Qt.ItemIsEditable)
            table.setItem(i, 0, key_item)
            
            # Value
            val_item = QTableWidgetItem(_summarize(val))
            val_item.setData(Qt.UserRole, val)
            val_item.setFlags(val_item.flags() & ~Qt.ItemIsEditable)
            table.setItem(i, 1, val_item)
            
            # Types
            types_item = QTableWidgetItem(f"K:{type(key).__name__}, V:{type(val).__name__}")
            types_item.setFlags(types_item.flags() & ~Qt.ItemIsEditable)
            table.setItem(i, 2, types_item)
    
    def on_cell_double_clicked(row, col):
        # Double-click to edit key or value
        active_entries = [e for i, e in enumerate(edited_entries) if i not in deleted_indices]
        if row >= len(active_entries):
            return
        
        key, val = active_entries[row]
        
        # Find the original index in edited_entries
        idx = 0
        active_idx = 0
        for i, (k, v) in enumerate(edited_entries):
            if i in deleted_indices:
                continue
            if active_idx == row:
                idx = i
                break
            active_idx += 1
        
        if col == 0:  # Edit key
            new_key = dispatch_edit(key, parent=dialog)
            if new_key is not None:
                edited_entries[idx] = (new_key, val)
                refresh_table()
        elif col == 1:  # Edit value
            new_val = dispatch_edit(val, parent=dialog)
            if new_val is not None:
                edited_entries[idx] = (key, new_val)
                refresh_table()
    
    def on_add_entry():
        """Add a new entry to the map"""
        # Try to infer key/value types from existing entries
        key_type = None
        val_type = None
        
        if len(edited_entries) > 0:
            first_entry = edited_entries[0]
            key_type = type(first_entry[0])
            val_type = type(first_entry[1])
        else:
            # For empty maps, try to infer from nanobind signatures
            inferred_keys, inferred_vals = _infer_map_types(value)
            if inferred_keys:
                key_type = _choose_type_dialog(inferred_keys, "Select key type", dialog)
            if inferred_vals:
                val_type = _choose_type_dialog(inferred_vals, "Select value type", dialog)
        
        # Create default instances
        try:
            if key_type is not None:
                new_key = _construct_default(key_type)
            else:
                QMessageBox.warning(dialog, "Cannot Add", "Cannot determine key type for this map")
                return
                
            if val_type is not None:
                new_val = _construct_default(val_type)
            else:
                QMessageBox.warning(dialog, "Cannot Add", "Cannot determine value type for this map")
                return
            if new_key is None or new_val is None:
                QMessageBox.warning(dialog, "Cannot Create", "Could not create default instances for the selected types")
                return
        except Exception as e:
            QMessageBox.warning(dialog, "Cannot Create", f"Could not create default instances: {e}")
            return
        
        # Edit the new key and value
        edited_key = dispatch_edit(new_key, parent=dialog)
        if edited_key is None:
            return
        
        edited_val = dispatch_edit(new_val, parent=dialog)
        if edited_val is None:
            return
        
        edited_entries.append((edited_key, edited_val))
        refresh_table()
    
    def on_remove_entry():
        """Remove the selected entry"""
        current_row = table.currentRow()
        if current_row < 0:
            QMessageBox.information(dialog, "No Selection", "Please select a row to remove")
            return
        
        # Find the original index
        active_idx = 0
        for i in range(len(edited_entries)):
            if i not in deleted_indices:
                if active_idx == current_row:
                    deleted_indices.add(i)
                    break
                active_idx += 1
        
        refresh_table()
    
    table.cellDoubleClicked.connect(on_cell_double_clicked)
    
    # Add/Remove buttons
    toolbar_layout = QHBoxLayout()
    add_btn = QPushButton("+ Add")
    add_btn.setToolTip("Add a new key-value pair")
    remove_btn = QPushButton("− Remove")
    remove_btn.setToolTip("Remove the currently selected entry")
    toolbar_layout.addWidget(add_btn)
    toolbar_layout.addWidget(remove_btn)
    toolbar_layout.addStretch()
    layout.addLayout(toolbar_layout)
    
    layout.addWidget(table)
    
    add_btn.clicked.connect(on_add_entry)
    remove_btn.clicked.connect(on_remove_entry)
    
    # Dialog buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    # Initial population
    refresh_table()
    
    if dialog.exec_() == QDialog.Accepted:
        # Build result map
        try:
            result = original_type()
            active_entries = [e for i, e in enumerate(edited_entries) if i not in deleted_indices]
            for key, val in active_entries:
                result[key] = val
            return result
        except Exception:
            return None
    return None


def create_atm_keys_table(keys_list, parent=None):
    """
    Create a table widget for editing atmospheric keys (used by AtmPoint and AtmField).
    
    This handles the variant key types: SpeciesEnum, SpeciesIsotope, AtmKey, 
    QuantumIdentifier, ScatteringSpeciesProperty.
    
    Parameters
    ----------
    keys_list : list of (key, value) tuples
        Initial list of (key, AtmData or Numeric) pairs
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    tuple of (QWidget, list, callable)
        Returns (widget, dict_items list, refresh_callback)
        - widget: The table widget to add to a layout
        - dict_items: Mutable list of (key, value) pairs
        - refresh_callback: Function to call to refresh the table display
    """
    from . import edit as dispatch_edit
    from pyarts3 import arts
    
    widget = QWidget(parent)
    layout = QVBoxLayout(widget)
    layout.setContentsMargins(0, 0, 0, 0)
    
    # Store as list of (key, value) pairs for editing
    dict_items = list(keys_list)
    
    # Table for atmospheric keys
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Key Type", "Key", "Value"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    table.setSelectionBehavior(QTableWidget.SelectRows)
    table.setSelectionMode(QTableWidget.SingleSelection)
    
    def refresh_table():
        table.setRowCount(len(dict_items))
        for row, (key, val) in enumerate(dict_items):
            # Key type
            key_type = type(key).__name__
            item = QTableWidgetItem(key_type)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 0, item)
            
            # Key
            item = QTableWidgetItem(str(key))
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 1, item)
            
            # Value - check if it's AtmData or just Numeric
            if hasattr(val, 'data_type'):  # AtmData
                value_str = f"AtmData ({val.data_type})"
            else:  # Numeric
                value_str = str(val)
            item = QTableWidgetItem(value_str)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 2, item)
    
    # Double-click to edit
    def on_cell_double_clicked(row, col):
        if row >= len(dict_items):
            return
        key, val = dict_items[row]
        
        if col == 2:  # Edit value
            result = dispatch_edit(val, parent=widget)
            if result is not None:
                dict_items[row] = (key, result)
                refresh_table()
        elif col == 1:  # Edit key
            result = dispatch_edit(key, parent=widget)
            if result is not None:
                dict_items[row] = (result, val)
                refresh_table()
    
    table.cellDoubleClicked.connect(on_cell_double_clicked)
    
    # Toolbar for add/remove
    toolbar = QToolBar()
    
    # Key type selector for adding new entries
    key_type_combo = QComboBox()
    key_type_combo.addItems([
        "AtmKey",
        "SpeciesEnum", 
        "SpeciesIsotope",
        "QuantumIdentifier",
        "ScatteringSpeciesProperty"
    ])
    toolbar.addWidget(QLabel("New key type:"))
    toolbar.addWidget(key_type_combo)
    
    # Add button
    add_action = QAction("+ Add", widget)
    
    def add_entry():
        key_type_name = key_type_combo.currentText()
        
        # Create a default key based on type
        try:
            if key_type_name == "AtmKey":
                new_key = arts.AtmKey.t  # temperature
            elif key_type_name == "SpeciesEnum":
                new_key = arts.SpeciesEnum.Water
            elif key_type_name == "SpeciesIsotope":
                new_key = arts.SpeciesIsotope("H2O-161")
            elif key_type_name == "QuantumIdentifier":
                new_key = arts.QuantumIdentifier("H2O-161")
            elif key_type_name == "ScatteringSpeciesProperty":
                QMessageBox.warning(widget, "Not Implemented", 
                                   "ScatteringSpeciesProperty add not yet implemented")
                return
            else:
                return
            
            # Edit the key before adding
            edited_key = dispatch_edit(new_key, parent=widget)
            if edited_key is not None:
                # Create default value (Numeric for AtmPoint, AtmData for AtmField)
                # Caller can specify what to use
                new_value = 0.0  # Default to Numeric
                dict_items.append((edited_key, new_value))
                refresh_table()
        except Exception as e:
            QMessageBox.warning(widget, "Error", f"Failed to create key: {e}")
    
    add_action.triggered.connect(add_entry)
    toolbar.addAction(add_action)
    
    # Remove button
    remove_action = QAction("− Remove", widget)
    remove_action.setToolTip("Remove selected entry")
    
    def remove_entry():
        current_row = table.currentRow()
        if current_row >= 0 and current_row < len(dict_items):
            dict_items.pop(current_row)
            refresh_table()
    
    remove_action.triggered.connect(remove_entry)
    toolbar.addAction(remove_action)
    
    layout.addWidget(toolbar)
    layout.addWidget(QLabel(f"Atmospheric keys. Double-click to edit."))
    layout.addWidget(table)
    
    return widget, dict_items, refresh_table


def create_subsurface_keys_table(keys_list, default_value_factory, parent=None):
    """
    Create a table widget for editing subsurface keys (used by SubsurfacePoint and SubsurfaceField).
    
    Handles key types: SubsurfaceKey, SubsurfacePropertyTag.
    
    Parameters
    ----------
    keys_list : list[(key, value)]
        Initial list of pairs (key, value)
    default_value_factory : callable
        A function returning a default value when adding a new entry (e.g., Numeric(0.0) or SubsurfaceData(Numeric(0.0)))
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    tuple(widget, dict_items, refresh_callback)
        widget: QWidget to add to layouts
        dict_items: mutable list of (key, value)
        refresh_callback: function to refresh table contents
    """
    from . import edit as dispatch_edit
    from pyarts3 import arts
    
    widget = QWidget(parent)
    layout = QVBoxLayout(widget)
    layout.setContentsMargins(0, 0, 0, 0)
    
    dict_items = list(keys_list)
    
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Key Type", "Key", "Value"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    table.setSelectionBehavior(QTableWidget.SelectRows)
    table.setSelectionMode(QTableWidget.SingleSelection)
    
    def refresh_table():
        table.setRowCount(len(dict_items))
        for row, (key, val) in enumerate(dict_items):
            # Key type
            key_type = type(key).__name__
            item = QTableWidgetItem(key_type)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 0, item)
            
            # Key
            item = QTableWidgetItem(str(key))
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 1, item)
            
            # Value
            item = QTableWidgetItem(str(val))
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 2, item)
    
    def on_cell_double_clicked(row, col):
        if row >= len(dict_items):
            return
        key, val = dict_items[row]
        
        if col == 2:  # Edit value
            result = dispatch_edit(val, parent=widget)
            if result is not None:
                dict_items[row] = (key, result)
                refresh_table()
        elif col == 1:  # Edit key
            result = dispatch_edit(key, parent=widget)
            if result is not None:
                dict_items[row] = (result, val)
                refresh_table()
    
    table.cellDoubleClicked.connect(on_cell_double_clicked)
    
    toolbar = QToolBar()
    key_type_combo = QComboBox()
    key_type_combo.addItems(["SubsurfaceKey", "SubsurfacePropertyTag"]) 
    toolbar.addWidget(QLabel("New key type:"))
    toolbar.addWidget(key_type_combo)
    
    add_action = QAction("+ Add", widget)
    
    def add_entry():
        key_type_name = key_type_combo.currentText()
        try:
            if key_type_name == "SubsurfaceKey":
                new_key = arts.SubsurfaceKey.temperature
            elif key_type_name == "SubsurfacePropertyTag":
                new_key = arts.SubsurfacePropertyTag("custom")
            else:
                return
            edited_key = dispatch_edit(new_key, parent=widget)
            if edited_key is not None:
                new_value = default_value_factory()
                dict_items.append((edited_key, new_value))
                refresh_table()
        except Exception as e:
            QMessageBox.warning(widget, "Error", f"Failed to create key: {e}")
    
    add_action.triggered.connect(add_entry)
    toolbar.addAction(add_action)
    
    remove_action = QAction("− Remove", widget)
    remove_action.setToolTip("Remove selected entry")
    
    def remove_entry():
        current_row = table.currentRow()
        if 0 <= current_row < len(dict_items):
            dict_items.pop(current_row)
            refresh_table()
    
    remove_action.triggered.connect(remove_entry)
    toolbar.addAction(remove_action)
    
    layout.addWidget(toolbar)
    layout.addWidget(QLabel("Subsurface keys. Double-click to edit."))
    layout.addWidget(table)
    
    return widget, dict_items, refresh_table


def create_surface_keys_table(keys_list, default_value_factory, parent=None):
    """
    Create a table widget for editing surface keys (used by SurfacePoint and SurfaceField).

    Handles key types: SurfaceKey, SurfacePropertyTag.

    Parameters
    ----------
    keys_list : list[(key, value)]
        Initial list of pairs (key, value)
    default_value_factory : callable
        A function returning a default value when adding a new entry (e.g., Numeric(0.0) or SurfaceData(Numeric(0.0)))
    parent : QWidget, optional
        Parent widget

    Returns
    -------
    tuple(widget, dict_items, refresh_callback)
        widget: QWidget to add to layouts
        dict_items: mutable list of (key, value)
        refresh_callback: function to refresh table contents
    """
    from . import edit as dispatch_edit
    from pyarts3 import arts

    widget = QWidget(parent)
    layout = QVBoxLayout(widget)
    layout.setContentsMargins(0, 0, 0, 0)

    dict_items = list(keys_list)

    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Key Type", "Key", "Value"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    table.setSelectionBehavior(QTableWidget.SelectRows)
    table.setSelectionMode(QTableWidget.SingleSelection)

    def refresh_table():
        table.setRowCount(len(dict_items))
        for row, (key, val) in enumerate(dict_items):
            # Key type
            key_type = type(key).__name__
            item = QTableWidgetItem(key_type)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 0, item)

            # Key
            item = QTableWidgetItem(str(key))
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 1, item)

            # Value
            item = QTableWidgetItem(str(val))
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 2, item)

    def on_cell_double_clicked(row, col):
        if row >= len(dict_items):
            return
        key, val = dict_items[row]

        if col == 2:  # Edit value
            result = dispatch_edit(val, parent=widget)
            if result is not None:
                dict_items[row] = (key, result)
                refresh_table()
        elif col == 1:  # Edit key
            result = dispatch_edit(key, parent=widget)
            if result is not None:
                dict_items[row] = (result, val)
                refresh_table()

    table.cellDoubleClicked.connect(on_cell_double_clicked)

    toolbar = QToolBar()
    key_type_combo = QComboBox()
    key_type_combo.addItems(["SurfaceKey", "SurfacePropertyTag"])  
    toolbar.addWidget(QLabel("New key type:"))
    toolbar.addWidget(key_type_combo)

    add_action = QAction("+ Add", widget)

    def add_entry():
        key_type_name = key_type_combo.currentText()
        try:
            if key_type_name == "SurfaceKey":
                # Choose a sensible default
                new_key = arts.SurfaceKey.t
            elif key_type_name == "SurfacePropertyTag":
                new_key = arts.SurfacePropertyTag("custom")
            else:
                return
            edited_key = dispatch_edit(new_key, parent=widget)
            if edited_key is not None:
                new_value = default_value_factory()
                dict_items.append((edited_key, new_value))
                refresh_table()
        except Exception as e:
            QMessageBox.warning(widget, "Error", f"Failed to create key: {e}")

    add_action.triggered.connect(add_entry)
    toolbar.addAction(add_action)

    remove_action = QAction("− Remove", widget)
    remove_action.setToolTip("Remove selected entry")

    def remove_entry():
        current_row = table.currentRow()
        if 0 <= current_row < len(dict_items):
            dict_items.pop(current_row)
            refresh_table()

    remove_action.triggered.connect(remove_entry)
    toolbar.addAction(remove_action)

    layout.addWidget(toolbar)
    layout.addWidget(QLabel("Surface keys. Double-click to edit."))
    layout.addWidget(table)

    return widget, dict_items, refresh_table
