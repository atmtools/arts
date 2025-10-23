"""Editor for XsecRecord type."""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem,
    QDialogButtonBox, QHeaderView, QLabel
)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit an XsecRecord object.
    
    XsecRecord contains cross-section data with fit coefficients:
    - species: SpeciesEnum
    - version: int
    - fitminpressures: Vector
    - fitmaxpressures: Vector
    - fitmintemperatures: Vector
    - fitmaxtemperatures: Vector
    - fitcoeffs: ArrayOfGriddedField1Named
    
    Parameters
    ----------
    value : XsecRecord
        The XsecRecord object to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    XsecRecord or None
        Modified XsecRecord if OK clicked, None if cancelled
    """
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit XsecRecord")
    dialog.setMinimumWidth(750)
    dialog.setMinimumHeight(500)
    
    layout = QVBoxLayout(dialog)
    
    layout.addWidget(QLabel("Cross-section absorption record. Double-click Value to edit."))
    
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Field", "Type", "Value"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    
    fields = [
        ("species", value.species, "SpeciesEnum"),
        ("version", value.version, "int"),
        ("fitminpressures", value.fitminpressures, "Vector"),
        ("fitmaxpressures", value.fitmaxpressures, "Vector"),
        ("fitmintemperatures", value.fitmintemperatures, "Vector"),
        ("fitmaxtemperatures", value.fitmaxtemperatures, "Vector"),
        ("fitcoeffs", value.fitcoeffs, "ArrayOfGriddedField1Named"),
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
            
            # Format value display
            if field_name == "fitcoeffs":
                value_str = f"[{len(val)} GriddedField1Named]"
            elif field_name in ["fitminpressures", "fitmaxpressures", "fitmintemperatures", "fitmaxtemperatures"]:
                value_str = f"Vector({len(val)} elements)"
            else:
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
        # XsecRecord: version is read-only static, other fields are writable
        result = arts.XsecRecord()
        for field_name, _, _ in fields:
            if field_name != "version":  # Skip read-only version field
                setattr(result, field_name, current_values[field_name])
        return result
    
    return None
