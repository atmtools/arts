"""Common reusable GUI components and utilities."""

import numpy as np
from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel, QPushButton,
    QDialogButtonBox, QSpinBox, QDoubleSpinBox, QLineEdit, QTextEdit, QWidget,
    QTableWidget, QTableWidgetItem
)
from PyQt5.QtCore import Qt


__all__ = [
    'create_simple_editor_dialog',
    'create_details_button',
    'create_description_dialog',
    'create_parameter_row_with_details',
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


def create_simple_editor_dialog(title, widget, value_getter, parent=None):
    """
    Create a simple dialog with a single input widget.
    
    Parameters
    ----------
    title : str
        Dialog window title
    widget : QWidget
        The input widget (QSpinBox, QLineEdit, etc.)
    value_getter : callable
        Function to call on widget to get the edited value
    parent : QWidget, optional
        Parent widget
    
    Returns
    -------
    value or None
        The edited value if accepted, None if cancelled
    """
    # Ensure a QApplication exists (avoid Qt crashes if none is running)
    from PyQt5.QtWidgets import QApplication
    app = QApplication.instance()
    if app is None:
        # Create an application with no argv to minimize side effects
        app = QApplication([])
    
    dialog = QDialog(parent)
    dialog.setWindowTitle(title)
    
    layout = QVBoxLayout()
    form = QFormLayout()
    
    form.addRow("Value:", widget)
    layout.addLayout(form)
    
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    if dialog.exec_() == QDialog.Accepted:
        return value_getter()
    return None


def create_description_dialog(title, name, type_info, description, parent=None):
    """
    Create a dialog showing parameter/variable description.
    
    Parameters
    ----------
    title : str
        Dialog window title
    name : str
        Parameter/variable name (displayed in bold)
    type_info : str or list
        Type information (single type or list of types)
    description : str
        Full description text
    parent : QWidget, optional
        Parent widget
    
    Returns
    -------
    QDialog
        The configured dialog (not yet shown)
    """
    dialog = QDialog(parent)
    dialog.setWindowTitle(title)
    dialog.setMinimumWidth(500)
    
    layout = QVBoxLayout()
    
    # Name
    name_label = QLabel(f"<b>{name}</b>")
    layout.addWidget(name_label)
    
    # Type(s)
    if isinstance(type_info, list):
        if len(type_info) == 1:
            type_label = QLabel(f"Type: <code>{type_info[0]}</code>")
        else:
            type_label = QLabel(f"Types: <code>{', '.join(type_info)}</code>")
    else:
        type_label = QLabel(f"Type: <code>{type_info}</code>")
    type_label.setWordWrap(True)
    layout.addWidget(type_label)
    
    # Description
    if len(description) > 500:
        # Use scrollable text edit for long descriptions
        desc_widget = QTextEdit()
        desc_widget.setReadOnly(True)
        desc_widget.setPlainText(description)
        desc_widget.setMinimumHeight(300)
        layout.addWidget(desc_widget)
    else:
        desc_label = QLabel(f"<br><b>Description:</b><br>{description}")
        desc_label.setWordWrap(True)
        desc_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(desc_label)
    
    # Close button
    close_btn = QDialogButtonBox(QDialogButtonBox.Ok)
    close_btn.accepted.connect(dialog.accept)
    layout.addWidget(close_btn)
    
    dialog.setLayout(layout)
    return dialog


def create_details_button(tooltip="Show details", max_width=30):
    """
    Create a standard "?" details button.
    
    Parameters
    ----------
    tooltip : str
        Tooltip text to display on hover
    max_width : int
        Maximum button width in pixels
    
    Returns
    -------
    QPushButton
        Configured details button
    """
    btn = QPushButton("?")
    btn.setMaximumWidth(max_width)
    btn.setToolTip(tooltip)
    return btn


def create_parameter_row_with_details(
    value_widget, 
    has_description=False,
    description_callback=None,
    additional_widgets=None
):
    """
    Create a parameter row with value widget and optional details button.
    
    Parameters
    ----------
    value_widget : QWidget
        The main widget displaying/editing the value
    has_description : bool
        Whether a description is available
    description_callback : callable, optional
        Function to call when details button is clicked
    additional_widgets : list of QWidget, optional
        Additional widgets to add to the row (e.g., edit button, type selector)
    
    Returns
    -------
    QWidget
        Row widget with horizontal layout containing all components
    """
    row_widget = QWidget()
    row_layout = QHBoxLayout()
    row_layout.setContentsMargins(0, 0, 0, 0)
    
    # Add any additional widgets first (e.g., type selector)
    if additional_widgets:
        for widget in additional_widgets:
            if widget is not None:
                row_layout.addWidget(widget)
    
    # Add the main value widget
    row_layout.addWidget(value_widget, stretch=1)
    
    # Add details button if description available
    if has_description and description_callback:
        details_btn = create_details_button()
        details_btn.clicked.connect(description_callback)
        row_layout.addWidget(details_btn)
    
    row_widget.setLayout(row_layout)
    return row_widget


def create_numeric_spinbox(value=0.0):
    """
    Create a configured QDoubleSpinBox for numeric editing.
    
    Parameters
    ----------
    value : float
        Initial value
    
    Returns
    -------
    QDoubleSpinBox
        Configured spinbox widget
    """
    spin_box = QDoubleSpinBox()
    spin_box.setRange(-1e308, 1e308)
    spin_box.setDecimals(10)
    spin_box.setValue(float(value))
    return spin_box


def edit_numeric(value, parent=None):
    """
    Edit a numeric value using a double spin box.
    
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
    spin_box = create_numeric_spinbox(value)
    
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
    Edit a list-like value (list, tuple, etc.) in a table.
    
    Parameters
    ----------
    value : list, tuple, or array-like
        The list-like value to edit
    parent : QWidget, optional
        Parent widget
    
    Returns
    -------
    list or None
        Edited value as a list if accepted, None if cancelled
    """
    # Ensure a QApplication exists (avoid Qt crashes if none is running)
    from PyQt5.QtWidgets import QApplication
    app = QApplication.instance()
    if app is None:
        # Create an application with no argv to minimize side effects
        app = QApplication([])
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit List")
    dialog.resize(600, 400)
    
    layout = QVBoxLayout()
    
    # Convert to list if needed
    if hasattr(value, '__iter__') and not isinstance(value, str):
        data = list(value)
    else:
        data = [value]
    
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
    # Ensure a QApplication exists (avoid Qt crashes if none is running)
    from PyQt5.QtWidgets import QApplication
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
    
    # Handle higher dimensions (3D+) - reshape to 2D with index labels
    else:
        dialog.resize(900, 650)
        
        from PyQt5.QtGui import QColor
        import itertools
        
        original_shape = data.shape
        
        # Handle empty arrays
        if data.size == 0:
            layout.addWidget(QLabel(f"Empty array with shape: {original_shape}"))
            layout.addWidget(QLabel("Cannot edit empty multi-dimensional arrays."))
            layout.addWidget(QLabel("Please reshape or populate the array first."))
            buttons = QDialogButtonBox(QDialogButtonBox.Ok)
            buttons.accepted.connect(dialog.accept)
            layout.addWidget(buttons)
            dialog.setLayout(layout)
            dialog.exec_()
            return None
        
        last_dim_size = original_shape[-1]  # Columns
        
        # Reshape: flatten all but last dimension
        # e.g., (3, 4, 5, 6) -> (3*4*5, 6)
        reshaped = data.reshape(-1, last_dim_size)
        num_rows = reshaped.shape[0]
        
        # Info label
        shape_str = ' × '.join(map(str, original_shape))
        info = QLabel(f"Shape: {shape_str} (reshaped to {num_rows} × {last_dim_size} for editing)")
        layout.addWidget(info)
        
        # Table
        table = QTableWidget()
        table.setRowCount(num_rows)
        table.setColumnCount(last_dim_size)
        
        # Set column headers (last dimension indices)
        table.setHorizontalHeaderLabels([str(i) for i in range(last_dim_size)])
        
        # Calculate multi-dimensional indices for row labels
        # For shape (3, 4, 5, 6), we need indices like [0,0,0], [0,0,1], ..., [0,0,4], [0,1,0], etc.
        outer_shape = original_shape[:-1]  # All dimensions except last
        row_labels = []
        
        # Generate all combinations of indices for outer dimensions
        for indices in itertools.product(*[range(s) for s in outer_shape]):
            label = '[' + ','.join(map(str, indices)) + ',:]'
            row_labels.append(label)
        
        table.setVerticalHeaderLabels(row_labels)
        
        # Track where each outer dimension index changes (for visual separators)
        separator_rows = set()
        if len(outer_shape) > 0:
            # Add separator after each change in first outer dimension
            stride = int(np.prod(outer_shape[1:])) if len(outer_shape) > 1 else 1
            for i in range(stride, num_rows, stride):
                separator_rows.add(i)
        
        # Fill data - store original values in UserRole
        for i in range(num_rows):
            for j in range(last_dim_size):
                item = QTableWidgetItem(f"{reshaped[i, j]:.6g}")
                item.setData(Qt.UserRole, float(reshaped[i, j]))  # Store original
                
                # Add visual separator for new outer index
                if i in separator_rows:
                    # Make the row slightly grayed to show boundary
                    item.setBackground(QColor(240, 240, 240))
                
                table.setItem(i, j, item)
        
        # Make vertical header wider to accommodate multi-index labels
        table.verticalHeader().setMinimumWidth(100)
        
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
            for i in range(num_rows):
                row = []
                for j in range(last_dim_size):
                    item = table.item(i, j)
                    row.append(item.data(Qt.UserRole))
                result.append(row)
            
            # Reshape back to original shape
            result_array = np.array(result).reshape(original_shape)
            # Reconstruct original type if needed
            if original_type != np.ndarray:
                try:
                    return original_type(result_array)
                except Exception:
                    return result_array
            return result_array
        return None


def edit_griddedfield(value, parent=None):
    """
    Edit a GriddedField* value.
    
    GriddedFields have:
    - grids: list of grid vectors (one per dimension)
    - gridnames: list of names for each grid dimension
    - dataname: name of the data field
    - data: the n-dimensional data array
    
    Parameters
    ----------
    value : GriddedField instance
        A gridded field object (GriddedField1, GriddedField2, etc.)
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    same-type-as-input or None
        The edited gridded field if accepted, None if cancelled
    """
    from PyQt5.QtWidgets import QApplication
    
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    
    original_type = type(value)
    
    # Remember if grids/gridnames were tuples
    grids_was_tuple = isinstance(value.grids, tuple) if hasattr(value, 'grids') else False
    gridnames_was_tuple = isinstance(value.gridnames, tuple) if hasattr(value, 'gridnames') else False
    
    dialog = QDialog(parent)
    dialog.setWindowTitle(f"Edit {original_type.__name__}")
    dialog.resize(800, 600)
    
    layout = QVBoxLayout()
    
    # Info label
    # Determine dimensions robustly
    dim = getattr(original_type, 'dim', None)
    if dim is None:
        try:
            dim = len(value.grids) if getattr(value, 'grids', None) is not None else 0
        except Exception:
            dim = 0
        # Fallback: parse trailing number from class name (e.g., GeodeticField3)
        if dim == 0:
            import re
            m = re.search(r"(\d+)$", original_type.__name__)
            if m:
                dim = int(m.group(1))
    if not dim:
        dim = 1
    
    # Prepare safe defaults if value is empty/uninitialized
    def _default_gridnames(n):
        if original_type.__name__ == 'GeodeticField3' and n == 3:
            return ["Altitude", "Latitude", "Longitude"]
        return [f"Dim {i+1}" for i in range(n)]
    
    def _ensure_list(v):
        try:
            return list(v) if v is not None else []
        except Exception:
            return []
    
    # Store edited grids and gridnames (convert to list if tuple), with fallbacks
    raw_grids = getattr(value, 'grids', None)
    edited_grids = _ensure_list(raw_grids)
    if len(edited_grids) == 0:
        edited_grids = [np.array([0.0])] * int(dim)
    raw_gridnames = getattr(value, 'gridnames', None)
    edited_gridnames = _ensure_list(raw_gridnames)
    if len(edited_gridnames) != len(edited_grids):
        edited_gridnames = _default_gridnames(len(edited_grids))
    
    # Compute a safe data array
    raw_data = getattr(value, 'data', None)
    try:
        data_shape = tuple(len(g) for g in edited_grids) if edited_grids else getattr(raw_data, 'shape', (0,))
    except Exception:
        data_shape = getattr(raw_data, 'shape', (0,))
    if not hasattr(raw_data, 'shape') or any(s == 0 for s in data_shape):
        safe_shape = tuple(max(1, int(len(g))) for g in edited_grids) or (1,)
        edited_data = [np.zeros(safe_shape)]
    else:
        edited_data = [raw_data]
    
    info = QLabel(f"GriddedField with {dim} dimension(s), data shape: {edited_data[0].shape}")
    info.setStyleSheet("color: #555; font-weight: bold;")
    layout.addWidget(info)
    
    # Data name
    dataname_layout = QHBoxLayout()
    dataname_layout.addWidget(QLabel("Data Name:"))
    dataname_edit = QLineEdit(value.dataname)
    dataname_layout.addWidget(dataname_edit)
    dataname_layout.addStretch()
    layout.addLayout(dataname_layout)
    
    # Grids section
    grids_label = QLabel("Grids (double-click 'Grid Name' or 'Values' to edit):")
    grids_label.setStyleSheet("font-weight: bold; margin-top: 10px;")
    layout.addWidget(grids_label)
    
    # edited_grids and edited_gridnames prepared above
    
    grids_table = QTableWidget()
    grids_table.setColumnCount(4)
    grids_table.setHorizontalHeaderLabels(["Dim", "Grid Name", "Size", "Values"])
    grids_table.horizontalHeader().setStretchLastSection(True)
    
    def refresh_grids_table():
        grids_table.setRowCount(len(edited_grids))
        for i, (grid, name) in enumerate(zip(edited_grids, edited_gridnames)):
            # Dimension
            dim_item = QTableWidgetItem(str(i))
            dim_item.setFlags(dim_item.flags() & ~Qt.ItemIsEditable)
            grids_table.setItem(i, 0, dim_item)
            
            # Grid name (editable)
            name_item = QTableWidgetItem(name)
            grids_table.setItem(i, 1, name_item)
            
            # Size
            size_item = QTableWidgetItem(str(len(grid)))
            size_item.setFlags(size_item.flags() & ~Qt.ItemIsEditable)
            grids_table.setItem(i, 2, size_item)
            
            # Values summary (click to edit)
            grid_array = np.array(grid)
            if len(grid_array) <= 5:
                summary = str(list(grid_array))
            else:
                summary = f"[{grid_array[0]}, {grid_array[1]}, ..., {grid_array[-1]}]"
            values_item = QTableWidgetItem(summary)
            values_item.setData(Qt.UserRole, grid)
            values_item.setFlags(values_item.flags() & ~Qt.ItemIsEditable)
            grids_table.setItem(i, 3, values_item)
    
    def on_grid_cell_changed(row, col):
        if col == 1:  # Grid name changed
            item = grids_table.item(row, 1)
            if item is not None:
                edited_gridnames[row] = item.text()
    
    def on_grid_double_clicked(row, col):
        if col == 3:  # Values column - edit the grid vector
            item = grids_table.item(row, 3)
            if item is None:
                return
            current_grid = item.data(Qt.UserRole)
            # Use edit_ndarraylike directly to avoid circular import issues
            new_grid = edit_ndarraylike(current_grid, parent=dialog)
            if new_grid is not None:
                edited_grids[row] = new_grid
                refresh_grids_table()
    
    grids_table.cellChanged.connect(on_grid_cell_changed)
    grids_table.cellDoubleClicked.connect(on_grid_double_clicked)
    
    refresh_grids_table()
    layout.addWidget(grids_table)
    
    # Data section
    data_label = QLabel("Data array (double-click to edit):")
    data_label.setStyleSheet("font-weight: bold; margin-top: 10px;")
    layout.addWidget(data_label)
    
    # edited_data prepared above; use list to allow mutation in nested function
    
    # Data info table (single row, double-click to edit)
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
        except Exception as e:
            print(f"Error creating gridded field: {e}")
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
    from PyQt5.QtWidgets import QApplication, QMessageBox
    
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
        original_idx = 0
        for i, (k, v) in enumerate(edited_entries):
            if i not in deleted_indices:
                if original_idx == row:
                    original_idx = i
                    break
                original_idx += 1
        
        if col == 0:  # Edit key
            new_key = edit_ndarraylike(key, parent=dialog)
            if new_key is not None:
                edited_entries[i] = (new_key, val)
                refresh_table()
        elif col == 1:  # Edit value
            new_val = edit_ndarraylike(val, parent=dialog)
            if new_val is not None:
                edited_entries[i] = (key, new_val)
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
        
        # Create default instances
        try:
            if key_type is not None:
                new_key = key_type()
            else:
                QMessageBox.warning(dialog, "Cannot Add", "Cannot determine key type for empty map")
                return
                
            if val_type is not None:
                new_val = val_type()
            else:
                QMessageBox.warning(dialog, "Cannot Add", "Cannot determine value type for empty map")
                return
        except Exception as e:
            QMessageBox.warning(dialog, "Cannot Create", f"Could not create default instances: {e}")
            return
        
        # Edit the new key and value
        edited_key = edit_ndarraylike(new_key, parent=dialog)
        if edited_key is None:
            return
        
        edited_val = edit_ndarraylike(new_val, parent=dialog)
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
        except Exception as e:
            print(f"Error creating map: {e}")
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
    from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem,
                                  QToolBar, QComboBox, QLabel, QMessageBox, QAction)
    from PyQt5.QtCore import Qt
    
    widget = QWidget(parent)
    layout = QVBoxLayout(widget)
    layout.setContentsMargins(0, 0, 0, 0)
    
    # Import arts here to avoid circular imports
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    
    # Store as list of (key, value) pairs for editing
    dict_items = list(keys_list)
    
    # Table for atmospheric keys
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Key Type", "Key", "Value"])
    from PyQt5.QtWidgets import QHeaderView
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
    from PyQt5.QtWidgets import (
        QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem,
        QToolBar, QComboBox, QLabel, QMessageBox, QAction, QHeaderView
    )
    from PyQt5.QtCore import Qt
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    
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
    from PyQt5.QtWidgets import (
        QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem,
        QToolBar, QComboBox, QLabel, QMessageBox, QAction, QHeaderView
    )
    from PyQt5.QtCore import Qt
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit

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
