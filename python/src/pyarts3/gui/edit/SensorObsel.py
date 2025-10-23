"""Editor for SensorObsel type."""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem,
    QDialogButtonBox, QHeaderView, QLabel
)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit a SensorObsel object.
    
    SensorObsel (Sensor Observation Element) has:
    - f_grid: AscendingGrid (frequency grid)
    - poslos: SensorPosLosVector (position and line-of-sight vectors)
    - weight_matrix: SparseStokvecMatrix (weight matrix for observations)
    
    Parameters
    ----------
    value : SensorObsel
        The SensorObsel object to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    SensorObsel or None
        Modified SensorObsel if OK clicked, None if cancelled
    """
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit SensorObsel")
    dialog.setMinimumWidth(700)
    
    layout = QVBoxLayout(dialog)
    
    layout.addWidget(QLabel("Sensor observation element properties. Double-click Value to edit."))
    
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Field", "Type", "Value"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    
    fields = [
        ("f_grid", value.f_grid, "AscendingGrid"),
        ("poslos", value.poslos, "SensorPosLosVector"),
        ("weight_matrix", value.weight_matrix, "StokvecMatrix"),
    ]
    
    table.setRowCount(len(fields))
    current_values = {}
    
    def refresh_table():
        for row, (field_name, field_value, field_type) in enumerate(fields):
            # Field name
            item = QTableWidgetItem(field_name)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 0, item)
            
            # Type
            item = QTableWidgetItem(field_type)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 1, item)
            
            # Value
            if field_name in current_values:
                val = current_values[field_name]
            else:
                val = field_value
                current_values[field_name] = val
            
            value_str = str(val)
            if len(value_str) > 100:
                value_str = value_str[:100] + "..."
            item = QTableWidgetItem(value_str)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 2, item)
    
    refresh_table()
    
    def on_cell_double_clicked(row, col):
        if col != 2:  # Only edit Value column
            return
        
        field_name = fields[row][0]
        current_val = current_values[field_name]
        
        result = dispatch_edit(current_val, parent=dialog)
        if result is not None:
            current_values[field_name] = result
            refresh_table()
    
    table.cellDoubleClicked.connect(on_cell_double_clicked)
    layout.addWidget(table)
    
    # OK/Cancel buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    if dialog.exec_() == QDialog.Accepted:
        # SensorObsel requires constructor with validated dimensions
        # weight_matrix must have shape (poslos.size, f_grid.size)
        try:
            result = arts.SensorObsel(
                current_values["f_grid"],
                current_values["poslos"],
                current_values["weight_matrix"]
            )
            return result
        except Exception as e:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.warning(
                dialog,
                "Invalid SensorObsel",
                f"Cannot create SensorObsel with these values:\n\n{str(e)}\n\n"
                "Note: weight_matrix must have shape (poslos_size, f_grid_size)"
            )
            return None
    
    return None
