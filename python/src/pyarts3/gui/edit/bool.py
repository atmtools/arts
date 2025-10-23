"""Editor for bool type using a checkbox."""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QCheckBox, QDialogButtonBox, QLabel
)


def edit(value, parent=None):
    """
    Edit a boolean value using a checkbox.
    
    Parameters
    ----------
    value : bool
        The current boolean value
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    bool or None
        The edited boolean value if OK clicked, None if cancelled
    """
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit Boolean")
    dialog.setMinimumWidth(300)
    
    layout = QVBoxLayout(dialog)
    
    # Label
    label = QLabel("Set the boolean value:")
    layout.addWidget(label)
    
    # Checkbox
    checkbox = QCheckBox("Enabled" if value else "Disabled")
    checkbox.setChecked(bool(value))
    
    def on_state_changed(state):
        checkbox.setText("Enabled" if state else "Disabled")
    
    checkbox.stateChanged.connect(on_state_changed)
    layout.addWidget(checkbox)
    
    # Buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    if dialog.exec_() == QDialog.Accepted:
        return checkbox.isChecked()
    
    return None
