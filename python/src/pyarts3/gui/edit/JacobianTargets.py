"""Editor for JacobianTargets type."""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem,
    QDialogButtonBox, QHeaderView, QLabel
)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit a JacobianTargets object.
    
    JacobianTargets contains arrays of different target types:
    - atm: ArrayOfAtmTargets (atmospheric targets)
    - surf: ArrayOfSurfaceTarget (surface targets)
    - subsurf: list (subsurface targets)
    - line: ArrayOfLineTarget (line targets)
    - sensor: ArrayOfSensorTarget (sensor targets)
    - error: ArrayOfErrorTarget (error targets)
    
    Parameters
    ----------
    value : JacobianTargets
        The JacobianTargets object to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    JacobianTargets or None
        Modified JacobianTargets if OK clicked, None if cancelled
    """
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit JacobianTargets")
    dialog.setMinimumWidth(700)
    
    layout = QVBoxLayout(dialog)
    
    layout.addWidget(QLabel("Jacobian targets. Double-click Value to edit arrays."))
    
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Field", "Type", "Value"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    
    fields = [
        ("atm", value.atm, "ArrayOfAtmTargets"),
        ("surf", value.surf, "ArrayOfSurfaceTarget"),
        ("subsurf", value.subsurf, "list"),
        ("line", value.line, "ArrayOfLineTarget"),
        ("sensor", value.sensor, "ArrayOfSensorTarget"),
        ("error", value.error, "ArrayOfErrorTarget"),
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
            
            # Value - show length for arrays
            if field_name in current_values:
                val = current_values[field_name]
            else:
                val = field_value
                current_values[field_name] = val
            
            try:
                value_str = f"[{len(val)} items]"
            except:
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
        result = arts.JacobianTargets()
        result.atm = current_values["atm"]
        result.surf = current_values["surf"]
        result.subsurf = current_values["subsurf"]
        result.line = current_values["line"]
        result.sensor = current_values["sensor"]
        result.error = current_values["error"]
        return result
    
    return None
