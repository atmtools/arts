"""Editor for QuantumIdentifier type.

QuantumIdentifier has:
- isot: SpeciesIsotope (the isotopologue)
- state: QuantumState (map of quantum numbers)
"""

from PyQt5.QtWidgets import (
    QApplication, QDialog, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem,
    QDialogButtonBox
)
from PyQt5.QtCore import Qt

__all__ = ["edit"]


def edit(value, parent=None):
    """
    Edit a QuantumIdentifier value.
    
    Parameters
    ----------
    value : QuantumIdentifier
        The quantum identifier to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    QuantumIdentifier or None
        The edited value if accepted, None if cancelled
    """
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    
    # Lazy import to avoid circular imports
    from . import edit as dispatch_edit
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit QuantumIdentifier")
    dialog.resize(600, 250)
    
    layout = QVBoxLayout()
    
    # Info
    info = QLabel(f"Quantum Identifier: {value}")
    info.setStyleSheet("color: #555; font-weight: bold;")
    layout.addWidget(info)
    
    # Store edited values
    edited_isot = [value.isot]
    edited_state = [value.state]
    
    # Table showing fields (double-click to edit)
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Field", "Type", "Value"])
    table.horizontalHeader().setStretchLastSection(True)
    table.setRowCount(2)
    
    def refresh_table():
        # Row 0: isot
        field_item = QTableWidgetItem("isot")
        field_item.setFlags(field_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(0, 0, field_item)
        
        type_item = QTableWidgetItem(type(edited_isot[0]).__name__)
        type_item.setFlags(type_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(0, 1, type_item)
        
        value_item = QTableWidgetItem(str(edited_isot[0]))
        value_item.setData(Qt.UserRole, edited_isot[0])
        value_item.setFlags(value_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(0, 2, value_item)
        
        # Row 1: state
        field_item = QTableWidgetItem("state")
        field_item.setFlags(field_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(1, 0, field_item)
        
        type_item = QTableWidgetItem(type(edited_state[0]).__name__)
        type_item.setFlags(type_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(1, 1, type_item)
        
        state_str = str(edited_state[0]) if edited_state[0] else "(empty)"
        value_item = QTableWidgetItem(state_str)
        value_item.setData(Qt.UserRole, edited_state[0])
        value_item.setFlags(value_item.flags() & ~Qt.ItemIsEditable)
        table.setItem(1, 2, value_item)
    
    def on_cell_double_clicked(row, col):
        if row == 0:  # isot
            new_isot = dispatch_edit(edited_isot[0], parent=dialog)
            if new_isot is not None:
                edited_isot[0] = new_isot
                refresh_table()
        elif row == 1:  # state
            new_state = dispatch_edit(edited_state[0], parent=dialog)
            if new_state is not None:
                edited_state[0] = new_state
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
            result = arts.QuantumIdentifier()
            result.isot = edited_isot[0]
            result.state = edited_state[0]
            return result
        except Exception as e:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.warning(dialog, "Error", f"Cannot create QuantumIdentifier: {e}")
            return None
    return None
