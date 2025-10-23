"""Editor for SubsurfaceField struct type."""

from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel,
                              QTableWidget, QTableWidgetItem, QDialogButtonBox,
                              QHeaderView, QToolBar, QComboBox, QAction, QMessageBox)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit a SubsurfaceField object.

    - bottom_depth: Numeric (depth of the bottom of the subsurface)
    - dict-like keys: SubsurfaceKey and SubsurfacePropertyTag -> SubsurfaceData
    """
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    from ..common import create_subsurface_keys_table

    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit SubsurfaceField")
    dialog.setMinimumWidth(800)
    dialog.setMinimumHeight(500)

    layout = QVBoxLayout(dialog)

    # Top section: bottom_depth
    top_layout = QHBoxLayout()
    top_layout.addWidget(QLabel("bottom_depth:"))
    current_bd = [value.bottom_depth]
    bd_label = QLabel(str(current_bd[0]))
    top_layout.addWidget(bd_label)

    from PyQt5.QtWidgets import QPushButton
    edit_bd_btn = QPushButton("Edit")

    def edit_bd():
        result = dispatch_edit(current_bd[0], parent=dialog)
        if result is not None:
            current_bd[0] = result
            bd_label.setText(str(result))

    edit_bd_btn.clicked.connect(edit_bd)
    top_layout.addWidget(edit_bd_btn)
    top_layout.addStretch()
    layout.addLayout(top_layout)

    # Keys entries
    keys_list = [(k, value[k]) for k in value.keys()]

    def _default_ssdata():
        return arts.SubsurfaceData(arts.Numeric(0.0))

    keys_widget, dict_items, refresh_keys = create_subsurface_keys_table(
        keys_list, default_value_factory=_default_ssdata, parent=dialog
    )
    refresh_keys()

    layout.addWidget(QLabel(f"Subsurface field entries ({len(dict_items)} total). Double-click to edit."))
    layout.addWidget(keys_widget)

    # OK/Cancel buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)

    if dialog.exec_() == QDialog.Accepted:
        # Construct new SubsurfaceField
        try:
            result = arts.SubsurfaceField(current_bd[0])
            for key, ss_data in dict_items:
                result[key] = ss_data
            return result
        except Exception as e:
            print(f"Error creating SubsurfaceField: {e}")
            return None

    return None
