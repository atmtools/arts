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
    dim = getattr(original_type, 'dim', len(value.grids))
    info = QLabel(f"GriddedField with {dim} dimension(s), data shape: {value.data.shape}")
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
    
    # Store edited grids and gridnames (convert to list if tuple)
    edited_grids = list(value.grids) if value.grids else []
    edited_gridnames = list(value.gridnames) if value.gridnames else []
    
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
    
    # Store edited data
    edited_data = [value.data]  # Use list to allow mutation in nested function
    
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
        summary = f"min={data_array.min():.3g}, max={data_array.max():.3g}, mean={data_array.mean():.3g}"
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
