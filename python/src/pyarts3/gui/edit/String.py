"""Editor for String values."""

from PyQt5.QtWidgets import QDialog, QVBoxLayout, QFormLayout, QLineEdit, QDialogButtonBox
from PyQt5.QtCore import Qt

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit a String value.
    
    Parameters
    ----------
    value : String or str
        The string value to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    str or None
        The edited value if accepted, None if cancelled
    """
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit String")
    
    layout = QVBoxLayout()
    form = QFormLayout()
    
    line_edit = QLineEdit()
    line_edit.setText(str(value))
    form.addRow("Value:", line_edit)
    
    layout.addLayout(form)
    
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    if dialog.exec_() == QDialog.Accepted:
        return line_edit.text()
    return None
