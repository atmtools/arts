"""Editor for SurfaceField struct type."""

from PyQt5.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QDialogButtonBox,
)


def edit(value, parent=None):
    """
    Edit a SurfaceField object.

    - ellipsoid: Vector2 (semi-major, semi-minor)
    - dict-like keys: SurfaceKey and SurfacePropertyTag -> SurfaceData
    """
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    from ..common import create_surface_keys_table

    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit SurfaceField")
    dialog.setMinimumWidth(800)
    dialog.setMinimumHeight(500)

    layout = QVBoxLayout(dialog)

    # Top section: ellipsoid
    top_layout = QHBoxLayout()
    top_layout.addWidget(QLabel("ellipsoid:"))
    current_ellipsoid = [value.ellipsoid]
    ell_label = QLabel(str(current_ellipsoid[0]))
    top_layout.addWidget(ell_label)

    edit_ell_btn = QPushButton("Edit")

    def edit_ell():
        result = dispatch_edit(current_ellipsoid[0], parent=dialog)
        if result is not None:
            current_ellipsoid[0] = result
            ell_label.setText(str(result))

    edit_ell_btn.clicked.connect(edit_ell)
    top_layout.addWidget(edit_ell_btn)
    top_layout.addStretch()
    layout.addLayout(top_layout)

    # Keys entries via unified keys() view
    keys_list = [(k, value[k]) for k in value.keys()]

    def _default_surfdata():
        return arts.SurfaceData(arts.Numeric(0.0))

    keys_widget, dict_items, refresh_keys = create_surface_keys_table(
        keys_list, default_value_factory=_default_surfdata, parent=dialog
    )
    refresh_keys()

    layout.addWidget(
        QLabel(
            f"Surface field entries ({len(dict_items)} total). Double-click to edit."
        )
    )
    layout.addWidget(keys_widget)

    # OK/Cancel buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)

    if dialog.exec_() == QDialog.Accepted:
        # Construct a new SurfaceField
        # Note: Python binding exposes __init__(planet: String). Use a neutral planet
        # then overwrite ellipsoid and keys. (Planet string is not persisted in the object.)
        try:
            result = arts.SurfaceField("Earth")
            result.ellipsoid = current_ellipsoid[0]
            for key, sd in dict_items:
                result[key] = sd
            return result
        except Exception as e:
            print(f"Error creating SurfaceField: {e}")
            return None

    return None
