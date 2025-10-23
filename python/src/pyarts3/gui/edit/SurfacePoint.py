"""Editor for SurfacePoint struct type."""

from PyQt5.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QLabel,
    QTableWidget,
    QTableWidgetItem,
    QGroupBox,
    QDialogButtonBox,
    QHeaderView,
)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit a SurfacePoint object with inline keys table.

    Core fields: temperature [K], elevation [m], normal [Vector2]
    Keys: SurfaceKey and SurfacePropertyTag -> Numeric values
    """
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    from ..common import create_surface_keys_table

    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit SurfacePoint")
    dialog.resize(800, 600)

    layout = QVBoxLayout(dialog)

    # Core fields section
    core_group = QGroupBox("Core Fields")
    core_layout = QVBoxLayout(core_group)

    core_table = QTableWidget()
    core_table.setColumnCount(2)
    core_table.setHorizontalHeaderLabels(["Field", "Value"])
    core_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
    core_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)

    core_fields = ["temperature", "elevation", "normal"]
    core_values = {f: getattr(value, f) for f in core_fields}

    def refresh_core_table():
        core_table.setRowCount(len(core_fields))
        for row, field in enumerate(core_fields):
            name_item = QTableWidgetItem(field)
            name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable)
            core_table.setItem(row, 0, name_item)

            val_item = QTableWidgetItem(str(core_values[field]))
            val_item.setFlags(val_item.flags() & ~Qt.ItemIsEditable)
            core_table.setItem(row, 1, val_item)

    def on_core_cell_double_clicked(row, col):
        if col != 1:
            return
        field = core_fields[row]
        result = dispatch_edit(core_values[field], parent=dialog)
        if result is not None:
            core_values[field] = result
            refresh_core_table()

    core_table.cellDoubleClicked.connect(on_core_cell_double_clicked)
    refresh_core_table()

    core_layout.addWidget(QLabel("Double-click value to edit."))
    core_layout.addWidget(core_table)
    layout.addWidget(core_group)

    # Keys section
    keys_group = QGroupBox("Surface Keys")
    keys_layout = QVBoxLayout(keys_group)

    # Collect existing key-value pairs using keys() + __getitem__
    keys_list = [(k, value[k]) for k in value.keys()]

    def _default_numeric():
        return arts.Numeric(0.0)

    keys_widget, dict_items, refresh_keys = create_surface_keys_table(
        keys_list, default_value_factory=_default_numeric, parent=dialog
    )
    refresh_keys()

    keys_layout.addWidget(keys_widget)
    layout.addWidget(keys_group)

    # Dialog buttons
    button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    button_box.accepted.connect(dialog.accept)
    button_box.rejected.connect(dialog.reject)
    layout.addWidget(button_box)

    if dialog.exec_():
        result = arts.SurfacePoint()
        for f in core_fields:
            setattr(result, f, core_values[f])
        for key, val in dict_items:
            result[key] = val
        return result

    return None
