"""Editor for QuantumValue type.

QuantumValue has:
- value: can be either String or Rational

Important: The type cannot be changed during editing.
"""

from PyQt5.QtWidgets import (
    QApplication, QDialog, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem,
    QDialogButtonBox
)
from PyQt5.QtCore import Qt

__all__ = ["edit"]


def edit(value, parent=None):
    """
    Edit a QuantumValue value.
    
    The type of the value (String or Rational) is preserved.
    
    Parameters
    ----------
    value : QuantumValue
        The quantum value to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    QuantumValue or None
        The edited value if accepted, None if cancelled
    """
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    
    # Lazy import to avoid circular imports
    from . import edit as dispatch_edit
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit QuantumValue")
    dialog.resize(500, 200)
    
    layout = QVBoxLayout()
    
    # Get the current value and its type
    current_value = value.value
    value_type = type(current_value).__name__
    
    # Info
    info = QLabel(f"Quantum Value: {value} (type: {value_type})")
    info.setStyleSheet("color: #555; font-weight: bold;")
    layout.addWidget(info)
    
    # Store edited value
    edited_value = [current_value]
    
    # Table showing the value (double-click to edit)
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Field", "Type", "Value"])
    table.horizontalHeader().setStretchLastSection(True)
    table.setRowCount(1)
    
    def refresh_table():
        # Field
        field_item = QTableWidgetItem("value")
        field_item.setFlags(field_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(0, 0, field_item)
        
        # Type
        type_item = QTableWidgetItem(type(edited_value[0]).__name__)
        type_item.setFlags(type_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(0, 1, type_item)
        
        # Value
        value_item = QTableWidgetItem(str(edited_value[0]))
        value_item.setData(Qt.UserRole, edited_value[0])
        value_item.setFlags(value_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(0, 2, value_item)
    
    def on_cell_double_clicked(row, col):
        new_value = dispatch_edit(edited_value[0], parent=dialog)
        if new_value is not None:
            # Check if type is preserved
            if type(new_value).__name__ != type(current_value).__name__:
                from PyQt5.QtWidgets import QMessageBox
                QMessageBox.warning(
                    dialog, 
                    "Type Mismatch", 
                    f"Type cannot be changed from {type(current_value).__name__} to {type(new_value).__name__}"
                )
                return
            edited_value[0] = new_value
            refresh_table()
    
    table.cellDoubleClicked.connect(on_cell_double_clicked)
    refresh_table()
    layout.addWidget(table)
    
    # Help text
    help_label = QLabel(f"Double-click to edit (type must remain {value_type})")
    help_label.setStyleSheet("color: #888; font-style: italic;")
    layout.addWidget(help_label)
    
    # Dialog buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    if dialog.exec_() == QDialog.Accepted:
        try:
            from pyarts3 import arts
            result = arts.QuantumValue(edited_value[0])
            return result
        except Exception as e:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.warning(dialog, "Error", f"Cannot create QuantumValue: {e}")
            return None
    return None
