"""Editor for Sun type."""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem,
    QDialogButtonBox, QHeaderView, QLabel
)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit a Sun object.
    
    Sun has the following attributes:
    - description: str (textual description)
    - distance: float (distance from Earth in meters)
    - latitude: float (latitude in degrees)
    - longitude: float (longitude in degrees)
    - radius: float (solar radius in meters)
    - spectrum: Matrix (solar spectrum data)
    
    Parameters
    ----------
    value : Sun
        The Sun object to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    Sun or None
        Modified Sun if OK clicked, None if cancelled
    """
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit Sun")
    dialog.setMinimumWidth(700)
    
    layout = QVBoxLayout(dialog)
    
    layout.addWidget(QLabel("Solar properties. Double-click Value to edit."))
    
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Field", "Type", "Value"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    
    fields = [
        ("description", value.description, "str"),
        ("distance", value.distance, "float"),
        ("latitude", value.latitude, "float"),
        ("longitude", value.longitude, "float"),
        ("radius", value.radius, "float"),
        ("spectrum", value.spectrum, "Matrix"),
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
        result = arts.Sun()
        result.description = current_values["description"]
        result.distance = current_values["distance"]
        result.latitude = current_values["latitude"]
        result.longitude = current_values["longitude"]
        result.radius = current_values["radius"]
        result.spectrum = current_values["spectrum"]
        return result
    
    return None
