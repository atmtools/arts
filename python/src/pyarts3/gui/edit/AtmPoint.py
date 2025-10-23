"""Editor for AtmPoint struct type."""

from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel, 
                              QTableWidget, QTableWidgetItem, QPushButton,
                              QGroupBox, QDialogButtonBox)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit an AtmPoint object with inline keys table.
    
    AtmPoint has core fields (temperature, pressure, mag, wind) and a
    keys() interface that returns atmospheric data values.
    
    Parameters
    ----------
    value : AtmPoint
        The AtmPoint object to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    AtmPoint or None
        Modified AtmPoint if OK clicked, None if cancelled
    """
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    from ..common import create_atm_keys_table
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit AtmPoint")
    dialog.resize(800, 600)
    
    layout = QVBoxLayout(dialog)
    
    # Core fields section
    core_group = QGroupBox("Core Fields")
    core_layout = QVBoxLayout(core_group)
    
    core_table = QTableWidget()
    core_table.setColumnCount(2)
    core_table.setHorizontalHeaderLabels(["Field", "Value"])
    from PyQt5.QtWidgets import QHeaderView
    core_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
    core_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
    
    # Core field names
    core_fields = ['temperature', 'pressure', 'mag', 'wind']
    core_values = {}
    
    for field in core_fields:
        core_values[field] = getattr(value, field)
    
    def refresh_core_table():
        core_table.setRowCount(len(core_fields))
        for row, field in enumerate(core_fields):
            # Field name
            name_item = QTableWidgetItem(field)
            name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable)
            core_table.setItem(row, 0, name_item)
            
            # Value
            val = core_values[field]
            val_item = QTableWidgetItem(str(val))
            val_item.setFlags(val_item.flags() & ~Qt.ItemIsEditable)
            core_table.setItem(row, 1, val_item)
    
    def on_core_cell_double_clicked(row, col):
        if col == 1:  # Edit value
            field = core_fields[row]
            old_val = core_values[field]
            result = dispatch_edit(old_val, parent=dialog)
            if result is not None:
                core_values[field] = result
                refresh_core_table()
    
    core_table.cellDoubleClicked.connect(on_core_cell_double_clicked)
    refresh_core_table()
    
    core_layout.addWidget(QLabel("Double-click value to edit."))
    core_layout.addWidget(core_table)
    layout.addWidget(core_group)
    
    # Keys section - use shared function
    keys_group = QGroupBox("Atmospheric Keys")
    keys_layout = QVBoxLayout(keys_group)
    
    # Get initial keys as list of (key, value) tuples
    keys_list = [(k, value[k]) for k in value.keys()]
    
    # Create the keys table widget
    keys_widget, dict_items, refresh_keys = create_atm_keys_table(keys_list, parent=dialog)
    refresh_keys()
    
    keys_layout.addWidget(keys_widget)
    layout.addWidget(keys_group)
    
    # Dialog buttons
    button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    button_box.accepted.connect(dialog.accept)
    button_box.rejected.connect(dialog.reject)
    layout.addWidget(button_box)
    
    if dialog.exec_():
        try:
            # Create new AtmPoint with updated core fields
            result = arts.AtmPoint()
            for field in core_fields:
                setattr(result, field, core_values[field])
            
            # Set all atmospheric keys
            for key, val in dict_items:
                result[key] = val
            
            return result
        except Exception as e:
            print(f"Error creating AtmPoint: {e}")
            return None
    
    return None
