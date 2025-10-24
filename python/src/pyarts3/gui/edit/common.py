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
