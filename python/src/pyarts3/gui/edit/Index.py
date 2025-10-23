"""Editor for Index (integer) values."""

from PyQt5.QtWidgets import QDialog, QVBoxLayout, QFormLayout, QSpinBox, QDialogButtonBox
from PyQt5.QtCore import Qt

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit an Index (integer) value.
    
    Parameters
    ----------
    value : Index or int
        The integer value to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    int or None
        The edited value if accepted, None if cancelled
    """
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit Index")
    
    layout = QVBoxLayout()
    form = QFormLayout()
    
    spin_box = QSpinBox()
    spin_box.setRange(-2147483648, 2147483647)
    spin_box.setValue(int(value))
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
