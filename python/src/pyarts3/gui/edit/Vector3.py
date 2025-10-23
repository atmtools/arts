"""Editor for Vector3 (3D vector) values."""

import numpy as np
from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QFormLayout, QDoubleSpinBox, 
                              QDialogButtonBox)
from PyQt5.QtCore import Qt

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit a Vector3 (3D vector) value.
    
    Parameters
    ----------
    value : Vector3 or 3-element array
        The 3D vector to edit (X, Y, Z components)
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    numpy.ndarray or None
        The edited value if accepted, None if cancelled
    """
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit Vector3")
    
    layout = QVBoxLayout()
    
    # Convert to numpy array if needed
    if hasattr(value, '__array__'):
        data = np.array(value)
    else:
        data = np.array(value)
    
    form = QFormLayout()
    spin_boxes = []
    
    labels = ["X", "Y", "Z"]
    
    for label, val in zip(labels, data):
        spin = QDoubleSpinBox()
        spin.setRange(-1e308, 1e308)
        spin.setDecimals(10)
        spin.setValue(float(val))
        form.addRow(f"{label}:", spin)
        spin_boxes.append(spin)
    
    layout.addLayout(form)
    
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    if dialog.exec_() == QDialog.Accepted:
        return np.array([spin.value() for spin in spin_boxes])
    return None
