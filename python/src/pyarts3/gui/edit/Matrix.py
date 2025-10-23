"""Editor for Matrix (2D array) values."""

import numpy as np
from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QLabel, QTableWidget, 
                              QTableWidgetItem, QDialogButtonBox, QScrollArea)
from PyQt5.QtCore import Qt

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit a Matrix (2D array) value.
    
    Parameters
    ----------
    value : Matrix or 2D array-like
        The matrix value to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    numpy.ndarray or None
        The edited value if accepted, None if cancelled
    """
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit Matrix")
    dialog.resize(800, 600)
    
    layout = QVBoxLayout()
    
    # Convert to numpy array if needed
    if hasattr(value, '__array__'):
        data = np.array(value)
    else:
        data = np.array(value)
    
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
    
    # Fill data
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            item = QTableWidgetItem(f"{data[i, j]:.6g}")
            table.setItem(i, j, item)
    
    # Make table scrollable
    scroll = QScrollArea()
    scroll.setWidget(table)
    scroll.setWidgetResizable(True)
    layout.addWidget(scroll)
    
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    if dialog.exec_() == QDialog.Accepted:
        # Build a nested Python list to feed back to original type constructor
        result = []
        for i in range(data.shape[0]):
            row = []
            for j in range(data.shape[1]):
                try:
                    row.append(float(table.item(i, j).text()))
                except (ValueError, AttributeError):
                    row.append(float(data[i, j]))  # Keep original on error
            result.append(row)
        # Preserve original ARTS type if possible
        try:
            return type(value)(result)
        except Exception:
            # Fallback to numpy array
            return np.array(result)
    return None
