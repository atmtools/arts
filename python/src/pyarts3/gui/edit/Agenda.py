"""Editor for ARTS Agenda values.

Provides a selector for workspace methods that can set the agenda,
plus an interface to call them with appropriate inputs.
"""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QDialogButtonBox, QApplication, QTextEdit, QPushButton,
    QMessageBox, QScrollArea, QWidget
)

__all__ = ["edit"]


def edit(value, parent=None):
    """
    Edit an ARTS Agenda by selecting and calling a workspace method.

    Parameters
    ----------
    value : Agenda
        An instance of an ARTS Agenda
    parent : QWidget, optional
        Parent widget for the dialog

    Returns
    -------
    Agenda or None
        The modified Agenda if accepted, None if cancelled
    """
    import pyarts3.arts as arts
    import pyarts3
    
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    
    # Get agenda name
    try:
        agenda_name = value.name
    except Exception:
        QMessageBox.warning(parent, "Error", "Cannot determine agenda name")
        return None
    
    # Get all workspace methods
    try:
        ws_methods = arts.globals.workspace_methods()
    except Exception as e:
        QMessageBox.warning(parent, "Error", f"Cannot retrieve workspace methods: {e}")
        return None
    
    # Filter methods that can set this agenda
    # A method can set the agenda if:
    # 1. Its output list is 1-long
    # 2. The output name matches agenda_name
    compatible_methods = []
    for method_name, method_info in ws_methods.items():
        try:
            outputs = method_info.output
            if len(outputs) == 1 and outputs[0] == agenda_name:
                compatible_methods.append(method_name)
        except Exception:
            continue
    
    if not compatible_methods:
        QMessageBox.information(
            parent,
            "No Compatible Methods",
            f"No workspace methods found that can set agenda '{agenda_name}'.\n\n"
            f"The agenda can only be viewed, not edited."
        )
        return None
    
    # Create dialog
    dialog = QDialog(parent)
    dialog.setWindowTitle(f"Edit Agenda: {agenda_name}")
    dialog.resize(700, 600)
    
    layout = QVBoxLayout()
    
    # Info section
    info_label = QLabel(f"<b>Agenda:</b> {agenda_name}")
    layout.addWidget(info_label)
    
    # Current agenda content (read-only preview)
    current_label = QLabel("<b>Current content:</b>")
    layout.addWidget(current_label)
    
    current_text = QTextEdit()
    current_text.setReadOnly(True)
    current_text.setMaximumHeight(150)
    try:
        current_text.setPlainText(repr(value))
    except Exception:
        current_text.setPlainText("<cannot display>")
    layout.addWidget(current_text)
    
    # Method selector
    method_label = QLabel(f"<b>Select method to set agenda:</b> ({len(compatible_methods)} available)")
    layout.addWidget(method_label)
    
    method_combo = QComboBox()
    method_combo.addItems(sorted(compatible_methods))
    layout.addWidget(method_combo)
    
    # Method description
    desc_label = QLabel("<b>Method description:</b>")
    layout.addWidget(desc_label)
    
    desc_text = QTextEdit()
    desc_text.setReadOnly(True)
    desc_text.setMaximumHeight(100)
    
    def update_description():
        selected = method_combo.currentText()
        if selected:
            try:
                method_info = ws_methods[selected]
                desc_parts = []
                if hasattr(method_info, 'desc'):
                    desc_parts.append(method_info.desc)
                if hasattr(method_info, 'input'):
                    inputs = method_info.input
                    if inputs:
                        desc_parts.append(f"\nInputs: {', '.join(inputs)}")
                desc_text.setPlainText('\n'.join(desc_parts) if desc_parts else "No description available")
            except Exception as e:
                desc_text.setPlainText(f"Error: {e}")
    
    method_combo.currentTextChanged.connect(update_description)
    update_description()
    layout.addWidget(desc_text)
    
    # Call method button
    call_button = QPushButton("Call Method to Set Agenda")
    
    def call_method():
        selected_method = method_combo.currentText()
        if not selected_method:
            return
        
        try:
            # Import the methods module to call workspace methods
            from . import methods
            
            # Create a temporary workspace
            ws = pyarts3.Workspace()
            
            # Get the method object
            method_obj = ws_methods[selected_method]
            
            # Call the method dialog to configure and execute
            methods.call_method_dialog(ws, selected_method, method_obj, parent=dialog)
            
            # Retrieve the modified agenda from workspace
            # Use eval as suggested: editted_agenda = ws.{agenda_name}
            modified_agenda = eval(f"ws.{agenda_name}")
            
            # Update preview
            try:
                current_text.setPlainText(repr(modified_agenda))
            except Exception:
                current_text.setPlainText("<cannot display>")
            
            # Store the modified agenda
            dialog.modified_agenda = modified_agenda
            
            QMessageBox.information(dialog, "Success", f"Method '{selected_method}' completed")
            
        except Exception as e:
            import traceback
            QMessageBox.critical(dialog, "Error", f"Failed to call method: {e}\n\n{traceback.format_exc()}")
    
    call_button.clicked.connect(call_method)
    layout.addWidget(call_button)
    
    # Store original agenda
    dialog.modified_agenda = value
    
    # Dialog buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    # Execute dialog
    if dialog.exec_() == QDialog.Accepted:
        return dialog.modified_agenda
    
    return None
