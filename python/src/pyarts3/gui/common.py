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
    
    # Convert to numpy array if needed
    if hasattr(value, '__array__'):
        data = np.array(value)
    else:
        data = np.array(value)
    
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
            return np.array(result)
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
            return np.array(result)
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
            return result_array
        return None
