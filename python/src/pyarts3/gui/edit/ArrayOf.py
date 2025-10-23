"""Generic editor for ArrayOf* (list-like) workspace types.

Displays a list of elements. Double-click an element to open its editor
(using the main edit dispatcher). Preserves original container type on save.
"""

from PyQt5.QtWidgets import (
    QApplication,
    QDialog, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem,
    QDialogButtonBox
)
from PyQt5.QtCore import Qt

__all__ = ["edit"]


def _summarize(value, max_len=100):
    try:
        s = str(value)
    except Exception:
        s = f"<{type(value).__name__}>"
    if len(s) > max_len:
        s = s[: max_len - 1] + "…"
    return s


def edit(value, parent=None):
    """
    Edit an ArrayOf* (list-like) value by editing individual elements.

    Parameters
    ----------
    value : ArrayOf* instance
        The list-like ARTS container to edit.
    parent : QWidget, optional
        Parent widget for the dialog

    Returns
    -------
    same-type-as-input or None
        The edited container if accepted, None if cancelled
    """
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])

    # Snapshot elements into a Python list for safe manipulation
    try:
        original = list(value)
    except Exception:
        # Fallback: try index-based access
        original = [x for x in value]

    n = len(original)
    edits = {}  # index -> new value

    # Lazy import to avoid circular imports
    from . import edit as dispatch_edit

    dialog = QDialog(parent)
    dialog.setWindowTitle(f"Edit {type(value).__name__}")
    dialog.resize(800, 500)

    layout = QVBoxLayout()

    info = QLabel(f"Size: {n}")
    info.setStyleSheet("color: #555;")
    layout.addWidget(info)

    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Index", "Type", "Value"])
    table.horizontalHeader().setStretchLastSection(True)
    table.setRowCount(n)

    for i, v in enumerate(original):
        # Index
        idx_it = QTableWidgetItem(str(i))
        idx_it.setFlags(idx_it.flags() & ~Qt.ItemIsEditable)
        table.setItem(i, 0, idx_it)
        # Type
        tname = type(v).__name__
        typ_it = QTableWidgetItem(tname)
        typ_it.setFlags(typ_it.flags() & ~Qt.ItemIsEditable)
        table.setItem(i, 1, typ_it)
        # Summary
        val_it = QTableWidgetItem(_summarize(v))
        val_it.setData(Qt.UserRole, v)  # store original object
        val_it.setFlags(val_it.flags() & ~Qt.ItemIsEditable)
        table.setItem(i, 2, val_it)

    def on_cell_double_clicked(row, col):
        # Allow double-click anywhere on the row to edit that element
        it = table.item(row, 2)
        if it is None:
            return
        current_obj = it.data(Qt.UserRole)
        original_type = type(current_obj)
        new_obj = dispatch_edit(current_obj, parent=dialog)
        if new_obj is not None:
            # Reconstruct original type if it changed
            result_type = type(new_obj)
            if original_type != result_type:
                try:
                    new_obj = original_type(new_obj)
                except Exception:
                    # If conversion fails, use new_obj as-is
                    pass
            edits[row] = new_obj
            # Update display
            it.setData(Qt.UserRole, new_obj)
            it.setText(_summarize(new_obj))
            table.item(row, 1).setText(type(new_obj).__name__)

    table.cellDoubleClicked.connect(on_cell_double_clicked)
    layout.addWidget(table)

    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)

    dialog.setLayout(layout)

    if dialog.exec_() == QDialog.Accepted:
        # Build final list, applying edits
        final = [edits.get(i, original[i]) for i in range(n)]
        try:
            return type(value)(final)
        except Exception:
            return final
    return None
