"""Generic read-only viewer for unsupported types."""

from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QLabel, QTextEdit, 
                              QDialogButtonBox)
from PyQt5.QtCore import Qt

__all__ = ['edit']


def edit(value, parent=None):
    """
    View (read-only) any value type.
    
    This is a fallback editor for types that don't have specific editors.
    It displays the value as text in a read-only view.
    
    Parameters
    ----------
    value : any
        The value to view
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    None
        Always returns None as this is read-only
    """
    dialog = QDialog(parent)
    dialog.setWindowTitle(f"View {type(value).__name__}")
    
    layout = QVBoxLayout()
    
    layout.addWidget(QLabel(f"Type: {type(value).__name__}"))
    layout.addWidget(QLabel("(Read-only view - no editor available)"))
    
    text_edit = QTextEdit()
    text_edit.setPlainText(str(value))
    text_edit.setReadOnly(True)
    layout.addWidget(text_edit)
    
    buttons = QDialogButtonBox(QDialogButtonBox.Ok)
    buttons.accepted.connect(dialog.accept)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    dialog.exec_()
    return None
