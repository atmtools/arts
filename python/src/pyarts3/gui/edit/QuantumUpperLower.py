"""Editor for QuantumUpperLower type.

QuantumUpperLower has:
- upper: QuantumValue
- lower: QuantumValue
"""

from PyQt5.QtWidgets import (
    QApplication, QDialog, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem,
    QDialogButtonBox
)
from PyQt5.QtCore import Qt

__all__ = ["edit"]


def edit(value, parent=None):
    """
    Edit a QuantumUpperLower value.
    
    Parameters
    ----------
    value : QuantumUpperLower
        The quantum upper/lower to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    QuantumUpperLower or None
        The edited value if accepted, None if cancelled
    """
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    
    # Lazy import to avoid circular imports
    from . import edit as dispatch_edit
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit QuantumUpperLower")
    dialog.resize(500, 200)
    
    layout = QVBoxLayout()
    
    # Info
    info = QLabel(f"Quantum Upper/Lower: {value}")
    info.setStyleSheet("color: #555; font-weight: bold;")
    layout.addWidget(info)
    
    # Store edited values
    edited_upper = [value.upper]
    edited_lower = [value.lower]
    
    # Table showing fields (double-click to edit)
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Field", "Type", "Value"])
    table.horizontalHeader().setStretchLastSection(True)
    table.setRowCount(2)
    
    def refresh_table():
        # Row 0: upper
        field_item = QTableWidgetItem("upper")
        field_item.setFlags(field_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(0, 0, field_item)
        
        type_item = QTableWidgetItem(type(edited_upper[0]).__name__)
        type_item.setFlags(type_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(0, 1, type_item)
        
        value_item = QTableWidgetItem(str(edited_upper[0]))
        value_item.setData(Qt.UserRole, edited_upper[0])
        value_item.setFlags(value_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(0, 2, value_item)
        
        # Row 1: lower
        field_item = QTableWidgetItem("lower")
        field_item.setFlags(field_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(1, 0, field_item)
        
        type_item = QTableWidgetItem(type(edited_lower[0]).__name__)
        type_item.setFlags(type_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(1, 1, type_item)
        
        value_item = QTableWidgetItem(str(edited_lower[0]))
        value_item.setData(Qt.UserRole, edited_lower[0])
        value_item.setFlags(value_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(1, 2, value_item)
    
    def on_cell_double_clicked(row, col):
        if row == 0:  # upper
            new_upper = dispatch_edit(edited_upper[0], parent=dialog)
            if new_upper is not None:
                edited_upper[0] = new_upper
                refresh_table()
        elif row == 1:  # lower
            new_lower = dispatch_edit(edited_lower[0], parent=dialog)
            if new_lower is not None:
                edited_lower[0] = new_lower
                refresh_table()
    
    table.cellDoubleClicked.connect(on_cell_double_clicked)
    refresh_table()
    layout.addWidget(table)
    
    # Help text
    help_label = QLabel("Double-click a field to edit it")
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
            result = arts.QuantumUpperLower()
            result.upper = edited_upper[0]
            result.lower = edited_lower[0]
            return result
        except Exception as e:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.warning(dialog, "Error", f"Cannot create QuantumUpperLower: {e}")
            return None
    return None
