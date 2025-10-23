"""Editor for CIARecord type."""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem,
    QDialogButtonBox, QHeaderView, QLabel
)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit a CIARecord object.
    
    CIARecord contains collision-induced absorption data:
    - specs: list (species pair specification)
    - data: ArrayOfGriddedField2 (CIA coefficient data)
    
    Parameters
    ----------
    value : CIARecord
        The CIARecord object to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    CIARecord or None
        Modified CIARecord if OK clicked, None if cancelled
    """
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit CIARecord")
    dialog.setMinimumWidth(650)
    dialog.setMinimumHeight(300)
    
    layout = QVBoxLayout(dialog)
    
    layout.addWidget(QLabel("Collision-induced absorption record. Double-click Value to edit."))
    
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Field", "Type", "Value"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    
    fields = [
        ("specs", value.specs, "list"),
        ("data", value.data, "ArrayOfGriddedField2"),
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
            
            if field_name == "specs":
                value_str = str(val)
            elif field_name == "data":
                value_str = f"[{len(val)} GriddedField2]"
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
        # CIARecord requires constructor: CIARecord(data, species1, species2)
        # specs is a tuple of two SpeciesEnum
        try:
            specs = current_values["specs"]
            data = current_values["data"]
            result = arts.CIARecord(data, specs[0], specs[1])
            return result
        except Exception as e:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.warning(dialog, "Error", f"Failed to create CIARecord: {e}")
            return None
    
    return None
