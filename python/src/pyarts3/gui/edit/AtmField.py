"""Editor for AtmField struct type."""

from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QTableWidget, QTableWidgetItem, 
                              QDialogButtonBox, QHeaderView, QLabel, QToolBar, QComboBox, 
                              QMessageBox, QPushButton, QAction)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit an AtmField object.
    
    AtmField has:
    - top_of_atmosphere: Numeric (float)
    - Atmospheric field dict with variant keys (SpeciesEnum, SpeciesIsotope, AtmKey, etc.)
    
    Parameters
    ----------
    value : AtmField
        The AtmField object to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    AtmField or None
        Modified AtmField if OK clicked, None if cancelled
    """
    from pyarts3 import arts
    
    # Lazy import to avoid circular imports
    from . import edit as dispatch_edit
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit AtmField")
    dialog.setMinimumWidth(800)
    dialog.setMinimumHeight(500)
    
    layout = QVBoxLayout(dialog)
    
    # Top section: top_of_atmosphere
    top_layout = QHBoxLayout()
    top_layout.addWidget(QLabel("top_of_atmosphere:"))
    current_toa = [value.top_of_atmosphere]
    toa_label = QLabel(str(current_toa[0]))
    top_layout.addWidget(toa_label)
    
    edit_toa_btn = QPushButton("Edit")
    
    def edit_toa():
        result = dispatch_edit(current_toa[0], parent=dialog)
        if result is not None:
            current_toa[0] = result
            toa_label.setText(str(result))
    
    edit_toa_btn.clicked.connect(edit_toa)
    top_layout.addWidget(edit_toa_btn)
    top_layout.addStretch()
    layout.addLayout(top_layout)
    
    # Get all keys from the AtmField
    all_keys = value.keys()
    
    # Store as list of (key, AtmData) pairs for editing
    # Access via value[key] to get the full AtmData object
    dict_items = [(key, value[key]) for key in all_keys]
    
    # Table for atmospheric field data
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Key Type", "Key", "Value (AtmData)"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    table.setSelectionBehavior(QTableWidget.SelectRows)
    table.setSelectionMode(QTableWidget.SingleSelection)
    
    def refresh_table():
        table.setRowCount(len(dict_items))
        for row, (key, val) in enumerate(dict_items):
            # Key type
            key_type = type(key).__name__
            item = QTableWidgetItem(key_type)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 0, item)
            
            # Key
            item = QTableWidgetItem(str(key))
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 1, item)
            
            # Value
            item = QTableWidgetItem(f"AtmData ({val.data_type})")
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 2, item)
    
    refresh_table()
    
    # Double-click to edit value
    def on_cell_double_clicked(row, col):
        if row >= len(dict_items):
            return
        key, val = dict_items[row]
        
        if col == 2:  # Edit value (AtmData)
            result = dispatch_edit(val, parent=dialog)
            if result is not None:
                dict_items[row] = (key, result)
                refresh_table()
        elif col == 1:  # Edit key
            result = dispatch_edit(key, parent=dialog)
            if result is not None:
                dict_items[row] = (result, val)
                refresh_table()
    
    table.cellDoubleClicked.connect(on_cell_double_clicked)
    
    # Toolbar for add/remove
    toolbar = QToolBar()
    
    # Key type selector for adding new entries
    key_type_combo = QComboBox()
    key_type_combo.addItems([
        "AtmKey",
        "SpeciesEnum", 
        "SpeciesIsotope",
        "QuantumIdentifier",
        "ScatteringSpeciesProperty"
    ])
    toolbar.addWidget(QLabel("New key type:"))
    toolbar.addWidget(key_type_combo)
    
    # Add button
    add_action = QAction("+ Add", dialog)
    
    def add_entry():
        key_type_name = key_type_combo.currentText()
        
        # Create a default key based on type
        try:
            if key_type_name == "AtmKey":
                new_key = arts.AtmKey.t  # temperature
            elif key_type_name == "SpeciesEnum":
                new_key = arts.SpeciesEnum.Water
            elif key_type_name == "SpeciesIsotope":
                new_key = arts.SpeciesIsotope("H2O-161")
            elif key_type_name == "QuantumIdentifier":
                new_key = arts.QuantumIdentifier("H2O-161")
            elif key_type_name == "ScatteringSpeciesProperty":
                # Need to check what's valid here
                QMessageBox.warning(dialog, "Not Implemented", 
                                   "ScatteringSpeciesProperty add not yet implemented")
                return
            else:
                return
            
            # Edit the key before adding
            edited_key = dispatch_edit(new_key, parent=dialog)
            if edited_key is not None:
                # Create default AtmData
                new_value = arts.AtmData()
                dict_items.append((edited_key, new_value))
                refresh_table()
        except Exception as e:
            QMessageBox.warning(dialog, "Error", f"Failed to create key: {e}")
    
    add_action.triggered.connect(add_entry)
    toolbar.addAction(add_action)
    
    # Remove button
    remove_action = QAction("− Remove", dialog)
    remove_action.setToolTip("Remove selected entry")
    
    def remove_entry():
        current_row = table.currentRow()
        if current_row >= 0 and current_row < len(dict_items):
            dict_items.pop(current_row)
            refresh_table()
    
    remove_action.triggered.connect(remove_entry)
    toolbar.addAction(remove_action)
    
    layout.addWidget(toolbar)
    layout.addWidget(QLabel(f"Atmospheric field entries ({len(dict_items)} total). Double-click to edit."))
    layout.addWidget(table)
    
    # OK/Cancel buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    if dialog.exec_() == QDialog.Accepted:
        # Create new AtmField with edited values
        result = arts.AtmField()
        result.top_of_atmosphere = current_toa[0]
        
        # Set each key-value pair
        for key, atm_data in dict_items:
            result[key] = atm_data
        
        return result
    
    return None
