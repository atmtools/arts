"""Editor for PropagationPathPoint type."""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem,
    QDialogButtonBox, QHeaderView, QLabel
)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit a PropagationPathPoint object.
    
    PropagationPathPoint has the following attributes:
    - los: Vector2 (line-of-sight angles: zenith, azimuth in degrees)
    - los_type: PathPositionType (position type for line-of-sight)
    - ngroup: float (group refractive index)
    - nreal: float (real refractive index)
    - pos: Vector3 (position: altitude, latitude, longitude)
    - pos_type: PathPositionType (position type for pos)
    
    Parameters
    ----------
    value : PropagationPathPoint
        The PropagationPathPoint object to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    PropagationPathPoint or None
        Modified PropagationPathPoint if OK clicked, None if cancelled
    """
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit PropagationPathPoint")
    dialog.setMinimumWidth(700)
    
    layout = QVBoxLayout(dialog)
    
    layout.addWidget(QLabel("Propagation path point properties. Double-click Value to edit."))
    
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Field", "Type", "Value"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    
    fields = [
        ("pos", value.pos, "Vector3"),
        ("pos_type", value.pos_type, "PathPositionType"),
        ("los", value.los, "Vector2"),
        ("los_type", value.los_type, "PathPositionType"),
        ("nreal", value.nreal, "float"),
        ("ngroup", value.ngroup, "float"),
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
        result = arts.PropagationPathPoint()
        result.pos = current_values["pos"]
        result.pos_type = current_values["pos_type"]
        result.los = current_values["los"]
        result.los_type = current_values["los_type"]
        result.nreal = current_values["nreal"]
        result.ngroup = current_values["ngroup"]
        return result
    
    return None
