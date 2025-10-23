"""Editor for Workspace values: browse and edit workspace variables."""

from PyQt5.QtWidgets import (
    QApplication,
    QDialog, QVBoxLayout, QHBoxLayout, QListWidget, QListWidgetItem,
    QLabel, QLineEdit, QDialogButtonBox, QMessageBox
)
from PyQt5.QtCore import Qt

__all__ = ['edit']


def edit(ws, parent=None):
    """
    Open a dialog listing workspace variables for browsing/editing.

    Parameters
    ----------
    ws : pyarts3.workspace.workspace.Workspace
        The workspace to browse and edit.
    parent : QWidget, optional
        Parent widget for the dialog.

    Returns
    -------
    Workspace or None
        Returns the workspace on close. No explicit edited value is returned.
    """
    # Ensure a QApplication exists (avoid Qt crashes if none is running)
    app = QApplication.instance()
    created_app = False
    if app is None:
        # Create an application with no argv to minimize side effects
        app = QApplication([])
        created_app = True

    # Lazy import to avoid circulars
    import pyarts3.gui.edit as editors
    import pyarts3.arts as cxx

    dialog = QDialog(parent)
    dialog.setWindowTitle("Workspace Variables")
    # Set minimum size to handle different DPI scaling across platforms
    dialog.setMinimumSize(800, 500)
    dialog.resize(900, 600)  # Wide enough for ~100 char variable names

    main = QVBoxLayout()

    # Header
    header = QLabel("Double-click a variable to view/edit its value")
    header.setStyleSheet("color: #555;")
    main.addWidget(header)

    # Search box
    search = QLineEdit()
    search.setPlaceholderText("Filter variables...")
    main.addWidget(search)

    # List of variables
    lst = QListWidget()
    lst.setSelectionMode(QListWidget.SingleSelection)
    main.addWidget(lst)

    # Footer info label
    info = QLabel("")
    info.setStyleSheet("color: #777;")
    main.addWidget(info)

    # Buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Close)
    buttons.rejected.connect(dialog.reject)  # Just in case
    buttons.accepted.connect(dialog.accept)
    buttons.button(QDialogButtonBox.Close).clicked.connect(dialog.accept)
    main.addWidget(buttons)

    dialog.setLayout(main)

    # Prepare variable metadata
    wsvars = cxx.globals.workspace_variables()  # name -> descriptor with .group
    all_wsv_names = list(wsvars.keys())  # Full list of all workspace variables

    def _group_name(name: str) -> str:
        """Return the declared workspace group/type name for variable name."""
        try:
            desc = wsvars.get(name)
            if not desc:
                return 'Unknown'
            g = getattr(desc, 'group', None)
            if g is None:
                g = getattr(desc, 'type', None)
            if isinstance(g, str):
                return g
            # If it's a type/class, return its __name__
            return getattr(g, '__name__', str(g))
        except Exception:
            return 'Unknown'

    # Populate list - only show variables that are set in ws.wsv
    def refresh_list(filter_text: str = ""):
        lst.clear()
        names = list(ws.wsv.keys())  # Only show set variables
        if filter_text:
            ft = filter_text.lower()
            names = [n for n in names if ft in n.lower()]
        names.sort()
        for name in names:
            # Determine type string and current status
            tname = _group_name(name)
            # Whether variable currently has a value
            has_value = False
            try:
                has_value = ws.has(name)
            except Exception:
                has_value = False
            status = "set" if has_value else "unset"
            text = f"{name}  —  {tname}  ({status})"
            item = QListWidgetItem(text)
            item.setData(Qt.UserRole, name)
            lst.addItem(item)
        # Add "+ Add workspace variable..." sentinel at bottom
        add_item = QListWidgetItem("+ Add workspace variable...")
        add_item.setData(Qt.UserRole, None)
        # Use a lighter gray that works across themes (Windows compatibility)
        from PyQt5.QtGui import QColor
        add_item.setForeground(QColor(128, 128, 128))
        lst.addItem(add_item)

        # Set info (exclude the sentinel from count)
        info.setText(f"{max(0, lst.count()-1)} variables")

    refresh_list()

    # Filtering behavior
    search.textChanged.connect(lambda s: refresh_list(s))

    # Add dialog to select and set an unset workspace variable
    def open_add_dialog():
        add_dlg = QDialog(dialog)
        add_dlg.setWindowTitle("Add workspace variable")
        add_dlg.setMinimumSize(800, 500)
        add_dlg.resize(900, 600)  # Match main dialog width for consistency
        add_layout = QVBoxLayout()

        add_info = QLabel("Double-click a variable to initialize and edit it")
        add_info.setStyleSheet("color: #555;")
        add_layout.addWidget(add_info)

        add_search = QLineEdit()
        add_search.setPlaceholderText("Filter variables...")
        add_layout.addWidget(add_search)

        add_list = QListWidget()
        add_layout.addWidget(add_list)

        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(add_dlg.reject)
        buttons.accepted.connect(add_dlg.accept)
        buttons.button(QDialogButtonBox.Close).clicked.connect(add_dlg.accept)
        add_layout.addWidget(buttons)

        add_dlg.setLayout(add_layout)

        def refresh_add_list(filter_text: str = ""):
            add_list.clear()
            # Show all workspace variables that are NOT in ws.wsv (not set)
            set_names = set(ws.wsv.keys())
            names = [n for n in all_wsv_names if n not in set_names]
            if filter_text:
                ft = filter_text.lower()
                names = [n for n in names if ft in n.lower()]
            names.sort()

            # Determine editor availability based on declared group name
            import pyarts3.gui.edit as editors_mod
            from PyQt5.QtGui import QColor
            for name in names:
                tname = _group_name(name)
                has_mod = hasattr(editors_mod, tname) and hasattr(getattr(editors_mod, tname), 'edit')
                text = f"{name}  —  {tname}"
                it = QListWidgetItem(text)
                it.setData(Qt.UserRole, (name, tname, has_mod))
                if not has_mod:
                    # Use explicit gray color for cross-platform consistency
                    it.setForeground(QColor(128, 128, 128))
                add_list.addItem(it)

        def on_add_double_clicked(it: QListWidgetItem):
            data = it.data(Qt.UserRole)
            if not data:
                return
            name, tname, has_mod = data
            if not has_mod:
                QMessageBox.information(add_dlg, "No editor",
                                        f"No editor available for type '{tname}'.")
                return
            # Initialize and edit
            try:
                ws.init(name)
                value = getattr(ws, name)
            except Exception as e:
                QMessageBox.warning(add_dlg, "Initialization failed",
                                    f"Could not initialize '{name}':\n{e}")
                return
            new_value = editors.edit(value, parent=add_dlg)
            if new_value is None:
                return
            try:
                setattr(ws, name, new_value)
            except Exception as e:
                QMessageBox.warning(add_dlg, "Set Failed",
                                    f"Failed to set '{name}':\n{e}")
                return
            # Close add dialog and refresh main list
            add_dlg.accept()
            refresh_list(search.text())

        add_list.itemDoubleClicked.connect(on_add_double_clicked)
        add_search.textChanged.connect(lambda s: refresh_add_list(s))

        refresh_add_list()
        add_dlg.exec_()

    # Double-click to edit selected variable or open add dialog
    def on_item_double_clicked(item: QListWidgetItem):
        name = item.data(Qt.UserRole)
        if not name:
            # Sentinel: open add dialog
            open_add_dialog()
            return
        # Resolve current value or create a default
        value = None
        try:
            if ws.has(name):
                value = getattr(ws, name)
            else:
                # Prefer initializing via the Workspace to avoid unsafe constructors
                try:
                    ws.init(name)
                    value = getattr(ws, name)
                except Exception:
                    # Fallback: try default-constructing the underlying type if available
                    tname = _group_name(name)
                    if tname and hasattr(cxx, tname):
                        try:
                            value = getattr(cxx, tname)()
                        except Exception:
                            QMessageBox.information(dialog, "Unavailable",
                                f"Cannot initialize variable '{name}' of type '{tname}'.")
                            return
                    else:
                        QMessageBox.information(dialog, "Unavailable",
                            f"Cannot determine type for '{name}'.")
                        return
        except Exception as e:
            QMessageBox.warning(dialog, "Error",
                                f"Failed to access '{name}':\n{e}")
            return

        # Open editor for the value
        new_value = editors.edit(value, parent=dialog)
        if new_value is None:
            return  # Cancelled or read-only

        # Assign back to workspace
        try:
            setattr(ws, name, new_value)
        except Exception as e:
            QMessageBox.warning(dialog, "Set Failed",
                                f"Failed to set '{name}':\n{e}")
            return

        # Update item display
        has_value = False
        try:
            has_value = ws.has(name)
        except Exception:
            has_value = True
        tname = _group_name(name)
        status = "set" if has_value else "unset"
        item.setText(f"{name}  —  {tname}  ({status})")

    lst.itemDoubleClicked.connect(on_item_double_clicked)

    dialog.exec_()

    # Do not quit the app if it existed before; if we created it, leave it running
    # (letting caller manage overall app lifecycle).  No explicit teardown here.
    return ws
