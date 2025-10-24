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


def edit_ndarraylike(value, parent=None, allowed_shape=None):
    """
    Edit an ndarray-like value (numpy array, ARTS arrays, etc.) in a table.
    Handles 1D and 2D arrays.
    
    Parameters
    ----------
    value : ndarray or array-like
        The array value to edit
    parent : QWidget, optional
        Parent widget
    allowed_shape : tuple[int, ...] | None, optional
        Shape constraint mask for resizing. Tuple length must equal array ndim.
        Use 0 for a modifiable dimension and a positive non-zero integer for a
        fixed dimension size. Example: (0, 3, 4) allows changing the first
        dimension only and enforces the last two to be 3 and 4.
        If None (default), shape cannot be changed in this editor.
    
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
    
    # Normalize and validate allowed_shape
    if allowed_shape is not None:
        try:
            allowed_shape = tuple(int(x) for x in allowed_shape)
        except Exception:
            allowed_shape = None
        if allowed_shape is not None and len(allowed_shape) != data.ndim:
            # Incompatible constraint; ignore but inform user
            QMessageBox.information(dialog, "Shape constraint ignored",
                                    "Provided allowed_shape length does not match array ndim; ignoring constraints.")
            allowed_shape = None

    # Handle 1D arrays
    if data.ndim == 1:
        dialog.resize(600, 400)
        
        # Info and shape controls
        info = QLabel("")
        def refresh_info():
            info.setText(f"Size: {len(data)}")
        refresh_info()
        layout.addWidget(info)

        # Optional shape control (length)
        length_spin = None
        if allowed_shape is not None:
            shape_row = QHBoxLayout()
            shape_row.addWidget(QLabel("Length:"))
            if allowed_shape[0] == 0:
                length_spin = QSpinBox()
                length_spin.setRange(1, 1_000_000)
                length_spin.setValue(int(len(data)) if len(data) > 0 else 1)
                shape_row.addWidget(length_spin)
            else:
                shape_row.addWidget(QLabel(str(int(allowed_shape[0]))))
            resize_btn = QPushButton("Resize")
            def do_resize_1d():
                nonlocal data
                new_len = len(data)
                if allowed_shape[0] == 0 and length_spin is not None:
                    new_len = int(length_spin.value())
                else:
                    # fixed
                    new_len = int(allowed_shape[0])
                if new_len <= 0:
                    QMessageBox.warning(dialog, "Invalid size", "Length must be >= 1")
                    return
                data = np.resize(data, (new_len,))
                refresh_info()
                build_table_1d()
            resize_btn.clicked.connect(do_resize_1d)
            shape_row.addWidget(resize_btn)
            shape_row.addStretch()
            layout.addLayout(shape_row)
        
        # Table
        table = QTableWidget()
        table.setColumnCount(2)
        table.setHorizontalHeaderLabels(["Index", "Value"])
        table.horizontalHeader().setStretchLastSection(True)

        def build_table_1d():
            table.setRowCount(len(data))
            for i, val in enumerate(data):
                # Index (read-only)
                idx_item = QTableWidgetItem(str(i))
                idx_item.setFlags(idx_item.flags() & ~Qt.ItemIsEditable)
                table.setItem(i, 0, idx_item)
                # Value (editable) - store original value in UserRole
                try:
                    v = float(val)
                except Exception:
                    v = np.nan
                val_item = QTableWidgetItem(f"{v:.10g}")
                val_item.setData(Qt.UserRole, float(v))  # Store original
                table.setItem(i, 1, val_item)
        build_table_1d()
        
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
        
        # Info label
        info = QLabel("")
        def refresh_info():
            info.setText(f"Shape: {data.shape[0]} × {data.shape[1]}")
        refresh_info()
        layout.addWidget(info)

        # Optional shape controls
        dim0_spin = None
        dim1_spin = None
        if allowed_shape is not None:
            shape_row = QHBoxLayout()
            shape_row.addWidget(QLabel("Shape:"))
            # Dim 0
            if allowed_shape[0] == 0:
                dim0_spin = QSpinBox()
                dim0_spin.setRange(1, 1_000_000)
                dim0_spin.setValue(int(data.shape[0]) if data.shape[0] > 0 else 1)
                shape_row.addWidget(dim0_spin)
            else:
                shape_row.addWidget(QLabel(str(int(allowed_shape[0]))))
            shape_row.addWidget(QLabel("×"))
            # Dim 1
            if allowed_shape[1] == 0:
                dim1_spin = QSpinBox()
                dim1_spin.setRange(1, 1_000_000)
                dim1_spin.setValue(int(data.shape[1]) if data.shape[1] > 0 else 1)
                shape_row.addWidget(dim1_spin)
            else:
                shape_row.addWidget(QLabel(str(int(allowed_shape[1]))))
            resize_btn = QPushButton("Resize")
            def do_resize_2d():
                nonlocal data
                new0 = data.shape[0]
                new1 = data.shape[1]
                new0 = int(dim0_spin.value()) if dim0_spin is not None else int(allowed_shape[0])
                new1 = int(dim1_spin.value()) if dim1_spin is not None else int(allowed_shape[1])
                if new0 <= 0 or new1 <= 0:
                    QMessageBox.warning(dialog, "Invalid size", "All dimensions must be >= 1")
                    return
                data = np.resize(data, (new0, new1))
                refresh_info()
                build_table_2d()
            resize_btn.clicked.connect(do_resize_2d)
            shape_row.addWidget(resize_btn)
            shape_row.addStretch()
            layout.addLayout(shape_row)
        
        # Table
        table = QTableWidget()
        def build_table_2d():
            # Reset table dimensions and headers
            table.setRowCount(int(data.shape[0]))
            table.setColumnCount(int(data.shape[1]))
            table.setHorizontalHeaderLabels([str(i) for i in range(int(data.shape[1]))])
            table.setVerticalHeaderLabels([str(i) for i in range(int(data.shape[0]))])
            # Fill data - store original values in UserRole
            for i in range(int(data.shape[0])):
                for j in range(int(data.shape[1])):
                    try:
                        v = float(data[i, j])
                    except Exception:
                        v = np.nan
                    item = QTableWidgetItem(f"{v:.6g}")
                    item.setData(Qt.UserRole, float(v))  # Store original
                    table.setItem(i, j, item)
        build_table_2d()
        
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
            for i in range(int(data.shape[0])):
                row = []
                for j in range(int(data.shape[1])):
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
        # Constrain data shape to current grid lengths
        try:
            grid_shape = tuple(len(np.array(g)) for g in edited_grids)
        except Exception:
            grid_shape = None
        # Use edit_ndarraylike directly with shape constraint
        new_data = edit_ndarraylike(current_data, parent=dialog, allowed_shape=grid_shape)
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


def edit_maplike(value, parent=None, allowed_key_types=None, allowed_value_types=None):
    """
    Edit a map-like value (dict-like with keys(), items(), values()).

        Contract (similar to edit_ndarraylike):
    - allowed_key_types is a sequence of classes/types for keys.
        - allowed_value_types is a sequence of classes/types for values.
    - If allowed_key_types is None: the editor is read-only. You can only view
      existing entries; no add, no remove, no editing.
        - If allowed_key_types is provided: you can add/remove entries and edit
      existing keys/values. When adding, you must select a key type from
            allowed_key_types; the value type is chosen as follows:
                1) If allowed_value_types is provided, select from that list
                2) Else, infer from the first existing entry
                3) If neither is available, adding is disabled

    Parameters
    ----------
    value : map-like instance
        Any dict-like object supporting items(), __setitem__, __delitem__
    parent : QWidget, optional
        Parent Qt widget
    allowed_key_types : Sequence[type] | None
        Allowed key types for adding new entries. None → read-only.
    allowed_value_types : Sequence[type] | None
        Allowed value types for adding new entries. If None, value type is
        inferred from the first existing entry.

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

    original_type = type(value)

    dialog = QDialog(parent)
    dialog.setWindowTitle(f"Edit {original_type.__name__}")
    dialog.resize(900, 600)

    layout = QVBoxLayout()

    # Info label
    info = QLabel(f"Map with {len(value)} entries")
    info.setStyleSheet("color: #555; font-weight: bold;")
    layout.addWidget(info)

    # Working copy as list of pairs for easier editing
    edited_entries = list(value.items()) if len(value) > 0 else []
    deleted_indices = set()

    # Table of entries
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
        return s if len(s) <= max_len else s[:max_len - 1] + "…"

    def refresh_table():
        active_entries = [e for i, e in enumerate(edited_entries) if i not in deleted_indices]
        info.setText(f"Map with {len(active_entries)} entries")
        table.setRowCount(len(active_entries))
        for i, (k, v) in enumerate(active_entries):
            ki = QTableWidgetItem(_summarize(k))
            ki.setData(Qt.UserRole, k)
            ki.setFlags(ki.flags() & ~Qt.ItemIsEditable)
            table.setItem(i, 0, ki)
            vi = QTableWidgetItem(_summarize(v))
            vi.setData(Qt.UserRole, v)
            vi.setFlags(vi.flags() & ~Qt.ItemIsEditable)
            table.setItem(i, 1, vi)
            ti = QTableWidgetItem(f"K:{type(k).__name__}, V:{type(v).__name__}")
            ti.setFlags(ti.flags() & ~Qt.ItemIsEditable)
            table.setItem(i, 2, ti)

    def _choose_type_dialog(type_list, title):
        if not type_list:
            return None
        if len(type_list) == 1:
            return type_list[0]
        dlg = QDialog(dialog)
        dlg.setWindowTitle(title)
        v = QVBoxLayout(dlg)
        v.addWidget(QLabel("Select a type:"))
        combo = QComboBox(dlg)
        combo.addItems([getattr(t, "__name__", str(t)) for t in type_list])
        v.addWidget(combo)
        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(dlg.accept)
        btns.rejected.connect(dlg.reject)
        v.addWidget(btns)
        return type_list[combo.currentIndex()] if dlg.exec_() == QDialog.Accepted else None

    def _construct_default(t):
        try:
            return t()
        except Exception:
            try:
                import pyarts3.arts as arts
                if t is getattr(arts, 'Numeric', None):
                    return t(0.0)
                if t is getattr(arts, 'Index', None):
                    return t(0)
                if t is getattr(arts, 'String', None):
                    return t("")
            except Exception:
                pass
            # Builtins
            if t in (int, float, str):
                return t()
        return None

    def on_cell_double_clicked(row, col):
        if allowed_key_types is None:
            return
        active_entries = [e for i, e in enumerate(edited_entries) if i not in deleted_indices]
        if row < 0 or row >= len(active_entries):
            return
        k, v = active_entries[row]
        # Find original index
        idx = 0
        active_idx = 0
        for i, (kk, vv) in enumerate(edited_entries):
            if i in deleted_indices:
                continue
            if active_idx == row:
                idx = i
                break
            active_idx += 1
        if col == 0:
            new_k = dispatch_edit(k, parent=dialog)
            if new_k is not None:
                edited_entries[idx] = (new_k, v)
                refresh_table()
        elif col == 1:
            new_v = dispatch_edit(v, parent=dialog)
            if new_v is not None:
                edited_entries[idx] = (k, new_v)
                refresh_table()

    def on_add_entry():
        if allowed_key_types is None:
            return
        # Choose key type
        key_types = list(allowed_key_types) if allowed_key_types is not None else []
        key_type = _choose_type_dialog(key_types, "Select key type")
        if key_type is None:
            return
        # Choose value type according to contract
        if allowed_value_types:
            val_type = _choose_type_dialog(list(allowed_value_types), "Select value type")
            if val_type is None:
                return
        else:
            if len(edited_entries) == 0:
                QMessageBox.warning(dialog, "Cannot Add", "Cannot determine value type (no allowed_value_types and no existing entries).")
                return
            val_type = type(edited_entries[0][1])
        new_key = _construct_default(key_type)
        new_val = _construct_default(val_type)
        if new_key is None or new_val is None:
            QMessageBox.warning(dialog, "Cannot Create", "Could not create default key/value instances.")
            return
        edited_key = dispatch_edit(new_key, parent=dialog)
        if edited_key is None:
            return
        edited_val = dispatch_edit(new_val, parent=dialog)
        if edited_val is None:
            return
        edited_entries.append((edited_key, edited_val))
        refresh_table()

    def on_remove_entry():
        if allowed_key_types is None:
            QMessageBox.information(dialog, "Read-only", "This map is read-only.")
            return
        current_row = table.currentRow()
        if current_row < 0:
            QMessageBox.information(dialog, "No Selection", "Please select a row to remove")
            return
        # Map current row to original index
        active_idx = 0
        for i in range(len(edited_entries)):
            if i not in deleted_indices:
                if active_idx == current_row:
                    deleted_indices.add(i)
                    break
                active_idx += 1
        refresh_table()

    # Hook up interactions
    if allowed_key_types is not None:
        table.cellDoubleClicked.connect(on_cell_double_clicked)

    # Toolbar
    toolbar = QHBoxLayout()
    if allowed_key_types is not None:
        add_btn = QPushButton("+ Add")
        add_btn.setToolTip("Add a new key-value pair")
        add_btn.clicked.connect(on_add_entry)
        toolbar.addWidget(add_btn)
        rm_btn = QPushButton("− Remove")
        rm_btn.setToolTip("Remove the selected entry")
        rm_btn.clicked.connect(on_remove_entry)
        toolbar.addWidget(rm_btn)
    toolbar.addStretch()
    layout.addLayout(toolbar)

    layout.addWidget(table)

    # Dialog buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)

    dialog.setLayout(layout)

    # Initial table
    refresh_table()

    if dialog.exec_() == QDialog.Accepted:
        try:
            result = original_type()
            active_entries = [e for i, e in enumerate(edited_entries) if i not in deleted_indices]
            for k, v in active_entries:
                result[k] = v
            return result
        except Exception:
            return None
    return None
