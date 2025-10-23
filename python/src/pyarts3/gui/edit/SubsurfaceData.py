"""Editor for SubsurfaceData struct type."""

from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem,
                              QDialogButtonBox, QHeaderView, QLabel, QComboBox, QPushButton, QMessageBox)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit a SubsurfaceData object.

    Same behavior as AtmData: data type and value editing, plus extrapolation settings.
    Double-click Type on the 'data' row to change type; double-click Value to edit.
    """
    from pyarts3 import arts
    from . import edit as dispatch_edit

    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit SubsurfaceData")
    dialog.setMinimumWidth(700)

    layout = QVBoxLayout(dialog)

    help_label = QLabel("Double-click Type (data row) to change type, or Value to edit it. Others: double-click Value.")
    help_label.setStyleSheet("color: gray; font-style: italic;")
    layout.addWidget(help_label)

    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Field", "Type", "Value"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

    current_data = [value.data]

    fields = [
        ("data", value.data, type(value.data).__name__, False),
        ("data_type", value.data_type, "String", True),
        ("alt_low", value.alt_low, "InterpolationExtrapolation", False),
        ("alt_upp", value.alt_upp, "InterpolationExtrapolation", False),
        ("lat_low", value.lat_low, "InterpolationExtrapolation", False),
        ("lat_upp", value.lat_upp, "InterpolationExtrapolation", False),
        ("lon_low", value.lon_low, "InterpolationExtrapolation", False),
        ("lon_upp", value.lon_upp, "InterpolationExtrapolation", False),
    ]

    table.setRowCount(len(fields))

    current_values = {}

    def refresh_table():
        for row, (field_name, field_value, field_type, readonly) in enumerate(fields):
            item = QTableWidgetItem(field_name)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            if readonly:
                item.setForeground(Qt.gray)
            table.setItem(row, 0, item)

            if field_name == "data":
                type_str = type(current_data[0]).__name__
            else:
                type_str = field_type
            item = QTableWidgetItem(type_str)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            if readonly:
                item.setForeground(Qt.gray)
            table.setItem(row, 1, item)

            if field_name == "data":
                val = current_data[0]
            elif field_name in current_values:
                val = current_values[field_name]
            else:
                val = field_value
                current_values[field_name] = val
            value_str = str(val)
            if len(value_str) > 100:
                value_str = value_str[:100] + "..."
            item = QTableWidgetItem(value_str)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            if readonly:
                item.setForeground(Qt.gray)
            table.setItem(row, 2, item)

    refresh_table()

    def on_cell_double_clicked(row, col):
        field_name, _, _, readonly = fields[row]

        if field_name == "data":
            if col == 1:
                # Type change
                type_dialog = QDialog(dialog)
                type_dialog.setWindowTitle("Change data type")
                from PyQt5.QtWidgets import QVBoxLayout as _QVBoxLayout
                type_layout = _QVBoxLayout(type_dialog)
                type_layout.addWidget(QLabel(f"Current type: {type(current_data[0]).__name__}"))
                type_combo = QComboBox()
                type_combo.addItems(["GeodeticField3", "Numeric", "NumericTernaryOperator"])
                type_combo.setCurrentText(type(current_data[0]).__name__)
                type_layout.addWidget(type_combo)
                buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
                type_layout.addWidget(buttons)
                buttons.accepted.connect(type_dialog.accept)
                buttons.rejected.connect(type_dialog.reject)
                if type_dialog.exec_() == QDialog.Accepted:
                    new_type = type_combo.currentText()
                    try:
                        if new_type == "Numeric":
                            new_data = arts.Numeric(0.0)
                        elif new_type == "GeodeticField3":
                            import numpy as np
                            new_data = arts.GeodeticField3(
                                name="Empty",
                                data=np.zeros((1, 1, 1)),
                                grid_names=["Altitude", "Latitude", "Longitude"],
                                grids=[np.array([0.0]), np.array([0.0]), np.array([0.0])],
                            )
                        elif new_type == "NumericTernaryOperator":
                            new_data = arts.NumericTernaryOperator()
                        else:
                            return
                        current_data[0] = new_data
                        refresh_table()
                    except Exception as e:
                        QMessageBox.warning(dialog, "Error", f"Failed to create {new_type}: {e}")
                return
            elif col == 2:
                # Edit current value
                result = dispatch_edit(current_data[0], parent=dialog)
                if result is not None:
                    current_data[0] = result
                    refresh_table()
                return
            else:
                return

        if readonly:
            QMessageBox.information(dialog, "Read-only Field", f"The '{field_name}' field is read-only and cannot be edited.")
            return
        if col != 2:
            return
        current_val = current_values[field_name]
        result = dispatch_edit(current_val, parent=dialog)
        if result is not None:
            current_values[field_name] = result
            refresh_table()

    table.cellDoubleClicked.connect(on_cell_double_clicked)
    layout.addWidget(table)

    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)

    if dialog.exec_() == QDialog.Accepted:
        result = arts.SubsurfaceData(current_data[0])
        result.alt_low = current_values["alt_low"]
        result.alt_upp = current_values["alt_upp"]
        result.lat_low = current_values["lat_low"]
        result.lat_upp = current_values["lat_upp"]
        result.lon_low = current_values["lon_low"]
        result.lon_upp = current_values["lon_upp"]
        return result

    return None
