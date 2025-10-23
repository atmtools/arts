"""Editor for AbsorptionLine type."""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem,
    QDialogButtonBox, QHeaderView, QLabel
)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit an AbsorptionLine object.
    
    AbsorptionLine contains spectroscopic parameters for a single line:
    - f0: float (center frequency)
    - a: float (line strength)
    - e0: float (lower state energy)
    - gu: float (upper state degeneracy)
    - gl: float (lower state degeneracy)
    - ls: LineShapeModel
    - z: ZeemanLineModel
    - qn: QuantumState
    
    Parameters
    ----------
    value : AbsorptionLine
        The AbsorptionLine object to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    AbsorptionLine or None
        Modified AbsorptionLine if OK clicked, None if cancelled
    """
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit AbsorptionLine")
    dialog.setMinimumWidth(700)
    dialog.setMinimumHeight(450)
    
    layout = QVBoxLayout(dialog)
    
    layout.addWidget(QLabel("Spectroscopic line parameters. Double-click Value to edit."))
    
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Field", "Type", "Value"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    
    fields = [
        ("f0", value.f0, "float"),
        ("a", value.a, "float"),
        ("e0", value.e0, "float"),
        ("gu", value.gu, "float"),
        ("gl", value.gl, "float"),
        ("ls", value.ls, "LineShapeModel"),
        ("z", value.z, "ZeemanLineModel"),
        ("qn", value.qn, "QuantumState"),
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
        result = arts.AbsorptionLine()
        for field_name, _, _ in fields:
            setattr(result, field_name, current_values[field_name])
        return result
    
    return None
