"""Editor for Vector (1D array) values."""

import numpy as np
from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QLabel, QTableWidget, 
                              QTableWidgetItem, QDialogButtonBox)
from PyQt5.QtCore import Qt

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit a Vector (1D array) value.
    
    Parameters
    ----------
    value : Vector or array-like
        The vector value to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    numpy.ndarray or None
        The edited value if accepted, None if cancelled
    """
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit Vector")
    dialog.resize(600, 400)
    
    layout = QVBoxLayout()
    
    # Convert to numpy array if needed
    if hasattr(value, '__array__'):
        data = np.array(value)
    else:
        data = np.array(value)
    
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
        
        # Value (editable)
        val_item = QTableWidgetItem(f"{val:.10g}")
        table.setItem(i, 1, val_item)
    
    layout.addWidget(table)
    
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    if dialog.exec_() == QDialog.Accepted:
        result = []
        for i in range(table.rowCount()):
            try:
                result.append(float(table.item(i, 1).text()))
            except (ValueError, AttributeError):
                result.append(float(data[i]))  # Keep original on error
        # Preserve original ARTS type if possible
        try:
            return type(value)(result)
        except Exception:
            # Fallback to numpy array
            return np.array(result)
    return None
