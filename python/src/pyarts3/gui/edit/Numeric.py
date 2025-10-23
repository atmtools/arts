"""Editor for Numeric (floating point) values."""

from PyQt5.QtWidgets import QDialog, QVBoxLayout, QFormLayout, QDoubleSpinBox, QDialogButtonBox
from PyQt5.QtCore import Qt

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit a Numeric (floating point) value.
    
    Parameters
    ----------
    value : Numeric or float
        The numeric value to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    float or None
        The edited value if accepted, None if cancelled
    """
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit Numeric")
    
    layout = QVBoxLayout()
    form = QFormLayout()
    
    spin_box = QDoubleSpinBox()
    spin_box.setRange(-1e308, 1e308)
    spin_box.setDecimals(10)
    spin_box.setValue(float(value))
    form.addRow("Value:", spin_box)
    
    layout.addLayout(form)
    
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    if dialog.exec_() == QDialog.Accepted:
        return spin_box.value()
    return None
