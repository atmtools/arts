"""Generic editor for ArrayOf* (list-like) workspace types.

Displays a list of elements. Double-click an element to open its editor
(using the main edit dispatcher). Preserves original container type on save.
"""

from PyQt5.QtWidgets import (
    QApplication,
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QTableWidget, QTableWidgetItem,
    QDialogButtonBox, QPushButton, QMessageBox
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

    # Track edits and deletions
    edits = {}  # index -> new value
    deleted = set()  # indices to delete
    added = []  # new elements to append

    # Track edits and deletions
    edits = {}  # index -> new value
    deleted = set()  # indices to delete
    added = []  # new elements to append

    # Lazy import to avoid circular imports
    from . import edit as dispatch_edit

    dialog = QDialog(parent)
    dialog.setWindowTitle(f"Edit {type(value).__name__}")
    dialog.resize(800, 500)

    layout = QVBoxLayout()

    def refresh_table():
        """Rebuild table from current state"""
        # Compute current working list
        working = [edits.get(i, original[i]) for i in range(len(original)) if i not in deleted]
        working.extend(added)
        
        n = len(working)
        info.setText(f"Size: {n}")
        
        table.setRowCount(n)
        for i, v in enumerate(working):
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

    info = QLabel()
    info.setStyleSheet("color: #555;")
    layout.addWidget(info)

    # Add/Remove buttons in a toolbar-style layout
    toolbar_layout = QHBoxLayout()
    add_btn = QPushButton("+ Add")
    add_btn.setToolTip("Add a new element to the list")
    remove_btn = QPushButton("− Remove")
    remove_btn.setToolTip("Remove the currently selected row")
    toolbar_layout.addWidget(add_btn)
    toolbar_layout.addWidget(remove_btn)
    toolbar_layout.addStretch()
    layout.addLayout(toolbar_layout)

    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Index", "Type", "Value"])
    table.horizontalHeader().setStretchLastSection(True)
    table.setSelectionBehavior(QTableWidget.SelectRows)
    table.setSelectionMode(QTableWidget.SingleSelection)

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
            
            # Determine if this is an added element or original
            working = [edits.get(i, original[i]) for i in range(len(original)) if i not in deleted]
            num_original_remaining = len(working)
            
            if row < num_original_remaining:
                # Find the actual original index
                original_idx = 0
                for i in range(len(original)):
                    if i not in deleted:
                        if original_idx == row:
                            edits[i] = new_obj
                            break
                        original_idx += 1
            else:
                # It's an added element
                added_idx = row - num_original_remaining
                added[added_idx] = new_obj
            
            refresh_table()

    def on_add_element():
        """Add a new element to the list"""
        # Try to determine the element type
        element_type = None
        
        # Check if we can infer element type from container type name
        container_type_name = type(value).__name__
        if container_type_name.startswith("ArrayOf"):
            element_type_name = container_type_name[7:]  # Remove "ArrayOf" prefix
            
            # Try to find the type in arts
            try:
                from pyarts3 import arts
                if hasattr(arts, element_type_name):
                    element_type = getattr(arts, element_type_name)
            except Exception:
                pass
        
        # If we have existing elements, use their type
        if element_type is None and len(original) > 0:
            element_type = type(original[0])
        
        if element_type is None and len(added) > 0:
            element_type = type(added[0])
        
        # Create a default instance
        if element_type is not None:
            try:
                new_element = element_type()
            except Exception:
                QMessageBox.warning(dialog, "Cannot Create Element", 
                                  f"Could not create a default instance of {element_type.__name__}")
                return
        else:
            QMessageBox.warning(dialog, "Cannot Determine Type", 
                              "Could not determine the element type for this container")
            return
        
        # Edit the new element
        edited_element = dispatch_edit(new_element, parent=dialog)
        if edited_element is not None:
            added.append(edited_element)
            refresh_table()

    def on_remove_element():
        """Remove the selected element"""
        current_row = table.currentRow()
        if current_row < 0:
            QMessageBox.information(dialog, "No Selection", "Please select a row to remove")
            return
        
        # Determine if this is an added element or original
        working = [edits.get(i, original[i]) for i in range(len(original)) if i not in deleted]
        num_original_remaining = len(working)
        
        if current_row < num_original_remaining:
            # Find the actual original index and mark as deleted
            original_idx = 0
            for i in range(len(original)):
                if i not in deleted:
                    if original_idx == current_row:
                        deleted.add(i)
                        break
                    original_idx += 1
        else:
            # It's an added element
            added_idx = current_row - num_original_remaining
            del added[added_idx]
        
        refresh_table()

    # Connect signals
    table.cellDoubleClicked.connect(on_cell_double_clicked)
    add_btn.clicked.connect(on_add_element)
    remove_btn.clicked.connect(on_remove_element)
    
    layout.addWidget(table)

    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)

    dialog.setLayout(layout)
    
    # Initial table population
    refresh_table()

    if dialog.exec_() == QDialog.Accepted:
        # Build final list: non-deleted originals (with edits) + added elements
        final = [edits.get(i, original[i]) for i in range(len(original)) if i not in deleted]
        final.extend(added)
        try:
            return type(value)(final)
        except Exception:
            return final
    return None
