"""GUI for calling ARTS workspace methods.

This module provides a graphical interface for browsing and calling ARTS workspace
methods. Methods are displayed with their input/output counts, and only methods
with all required workspace variable inputs available are enabled for calling.

The method call dialog shows:
- Workspace Inputs: Read-only display of workspace variables used as inputs
- Generic Inputs (GIN): Editable parameters with default values
- Workspace Outputs: Variables that will be modified or created
- Generic Outputs (GOUT): New variables that will be created and need names

Usage:
    from pyarts3.gui.methods import show_methods_dialog
    from pyarts3.workspace import Workspace
    
    ws = Workspace()
    # ... set up workspace variables ...
    show_methods_dialog(ws)
"""

from PyQt5.QtWidgets import (
    QApplication, QDialog, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit,
    QListWidget, QListWidgetItem, QDialogButtonBox, QMessageBox, QWidget,
    QScrollArea, QFormLayout, QGroupBox, QPushButton
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor


__all__ = ['show_methods_dialog']


def show_methods_dialog(ws, parent=None):
    """
    Open a dialog to browse and call workspace methods.
    
    Parameters
    ----------
    ws : pyarts3.workspace.workspace.Workspace
        The workspace to operate on.
    parent : QWidget, optional
        Parent widget for the dialog.
    """
    app = QApplication.instance()
    if app is None:
        app = QApplication([])

    from pyarts3.arts.globals import workspace_methods

    dialog = QDialog(parent)
    dialog.setWindowTitle("Workspace Methods")
    dialog.setMinimumSize(900, 600)
    dialog.resize(1000, 700)

    main = QVBoxLayout()

    # Header
    header = QLabel("Select a method to call")
    header.setStyleSheet("color: #555;")
    main.addWidget(header)

    # Controls row: search box and checkbox
    controls_layout = QHBoxLayout()
    
    search = QLineEdit()
    search.setPlaceholderText("Filter methods...")
    controls_layout.addWidget(search)
    
    from PyQt5.QtWidgets import QCheckBox
    show_all_checkbox = QCheckBox("Show uncallable methods")
    show_all_checkbox.setChecked(False)
    controls_layout.addWidget(show_all_checkbox)
    
    main.addLayout(controls_layout)

    # List of methods
    methods_list = QListWidget()
    methods_list.setSelectionMode(QListWidget.SingleSelection)
    main.addWidget(methods_list)

    # Info label
    info = QLabel("")
    info.setStyleSheet("color: #777;")
    main.addWidget(info)

    # Buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Close)
    buttons.button(QDialogButtonBox.Close).clicked.connect(dialog.accept)
    main.addWidget(buttons)

    dialog.setLayout(main)

    # Get all methods
    all_methods = workspace_methods()
    # Sort case-insensitively for better alphabetical ordering
    method_names = sorted(all_methods.keys(), key=str.lower)

    def can_call_method(method_obj):
        """Check if all required inputs are available in the workspace."""
        # Check workspace variable inputs
        for inp in method_obj.input:
            if not ws.has(inp):
                return False
        # GIN inputs with no default value are required
        for i, gin_val in enumerate(method_obj.gin_value):
            if gin_val is None:
                # No default, must be provided by user (always allow for now)
                pass
        return True

    def refresh_methods(filter_text: str = ""):
        methods_list.clear()
        names = method_names
        if filter_text:
            ft = filter_text.lower()
            names = [n for n in names if ft in n.lower()]

        show_uncallable = show_all_checkbox.isChecked()
        enabled_count = 0
        total_count = 0
        
        for name in names:
            method_obj = all_methods[name]
            can_call = can_call_method(method_obj)
            
            # Skip uncallable methods unless checkbox is checked
            if not can_call and not show_uncallable:
                continue
            
            total_count += 1
            
            # Show input/output counts
            n_in = len(method_obj.input)
            n_out = len(method_obj.output)
            n_gin = len(method_obj.gin)
            n_gout = len(method_obj.gout)
            
            text = f"{name}  (in:{n_in}, out:{n_out}, gin:{n_gin}, gout:{n_gout})"
            item = QListWidgetItem(text)
            item.setData(Qt.UserRole, name)
            
            if not can_call:
                # Gray out uncallable methods when shown
                item.setForeground(QColor(128, 128, 128))
                item.setFlags(item.flags() & ~Qt.ItemIsEnabled)
            else:
                enabled_count += 1
                
            methods_list.addItem(item)

        if show_uncallable:
            info.setText(f"{total_count} methods ({enabled_count} callable)")
        else:
            info.setText(f"{enabled_count} callable methods")

    refresh_methods()
    search.textChanged.connect(lambda s: refresh_methods(s))
    show_all_checkbox.stateChanged.connect(lambda: refresh_methods(search.text()))

    def on_method_double_clicked(item: QListWidgetItem):
        method_name = item.data(Qt.UserRole)
        if not method_name:
            return
        
        method_obj = all_methods[method_name]
        if not can_call_method(method_obj):
            QMessageBox.warning(dialog, "Cannot call",
                              f"Method '{method_name}' cannot be called: missing required inputs")
            return
        
        # Open method call dialog
        call_method_dialog(ws, method_name, method_obj, dialog)
        # Refresh list in case workspace state changed
        refresh_methods(search.text())

    methods_list.itemDoubleClicked.connect(on_method_double_clicked)

    dialog.exec_()


def call_method_dialog(ws, method_name, method_obj, parent=None):
    """
    Dialog to configure and call a workspace method.
    
    Parameters
    ----------
    ws : Workspace
        The workspace instance.
    method_name : str
        Name of the method.
    method_obj : WorkspaceMethodInternalRecord
        The method record.
    parent : QWidget, optional
        Parent widget.
    """
    import pyarts3.gui.edit as editors
    import pyarts3.arts as cxx
    
    dialog = QDialog(parent)
    dialog.setWindowTitle(f"Call: {method_name}")
    dialog.setMinimumSize(800, 600)

    main = QVBoxLayout()

    # Method description header with details button
    desc_header = QWidget()
    desc_header_layout = QHBoxLayout()
    desc_header_layout.setContentsMargins(0, 0, 0, 0)
    
    desc_label = QLabel(f"<b>{method_name}</b>")
    desc_header_layout.addWidget(desc_label, stretch=1)
    
    # Full description button
    if len(method_obj.desc) > 200:
        full_desc_btn = QPushButton("Full Description...")
        full_desc_btn.setMaximumWidth(150)
        
        def show_full_desc():
            desc_dialog = QDialog(dialog)
            desc_dialog.setWindowTitle(f"Method: {method_name}")
            desc_dialog.setMinimumSize(700, 500)
            
            desc_dialog_layout = QVBoxLayout()
            
            method_label = QLabel(f"<h2>{method_name}</h2>")
            desc_dialog_layout.addWidget(method_label)
            
            # Use a scrollable text area for long descriptions
            from PyQt5.QtWidgets import QTextEdit
            desc_text_edit = QTextEdit()
            desc_text_edit.setReadOnly(True)
            desc_text_edit.setPlainText(method_obj.desc)
            desc_dialog_layout.addWidget(desc_text_edit)
            
            close_btn = QDialogButtonBox(QDialogButtonBox.Ok)
            close_btn.accepted.connect(desc_dialog.accept)
            desc_dialog_layout.addWidget(close_btn)
            
            desc_dialog.setLayout(desc_dialog_layout)
            desc_dialog.exec_()
        
        full_desc_btn.clicked.connect(show_full_desc)
        desc_header_layout.addWidget(full_desc_btn)
    
    desc_header.setLayout(desc_header_layout)
    main.addWidget(desc_header)
    
    desc_text = QLabel(method_obj.desc[:200] + ("..." if len(method_obj.desc) > 200 else ""))
    desc_text.setWordWrap(True)
    desc_text.setStyleSheet("color: #555;")
    main.addWidget(desc_text)

    # Scrollable area for inputs/outputs
    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll_widget = QWidget()
    scroll_layout = QVBoxLayout()
    scroll_widget.setLayout(scroll_layout)
    scroll.setWidget(scroll_widget)
    main.addWidget(scroll)

    # --- Workspace Inputs (read-only display) ---
    if method_obj.input:
        # Get workspace variable descriptions
        wsvars = cxx.globals.workspace_variables()
        
        input_group = QGroupBox("Workspace Inputs (from workspace)")
        input_layout = QFormLayout()
        for inp in method_obj.input:
            # Row with value display and details button
            input_row = QWidget()
            input_row_layout = QHBoxLayout()
            input_row_layout.setContentsMargins(0, 0, 0, 0)
            
            value_str = "<not set>"
            if ws.has(inp):
                try:
                    val = getattr(ws, inp)
                    value_str = f"<{type(val).__name__}>"
                except Exception as e:
                    value_str = f"<error: {e}>"
            
            label = QLabel(value_str)
            label.setStyleSheet("color: #555;")
            input_row_layout.addWidget(label, stretch=1)
            
            # Add "?" button if variable has description
            if inp in wsvars:
                wsvar = wsvars[inp]
                if hasattr(wsvar, 'desc') and wsvar.desc:
                    input_details_btn = QPushButton("?")
                    input_details_btn.setMaximumWidth(30)
                    input_details_btn.setToolTip("Show variable description")
                    
                    def make_input_details_handler(var_name, var_obj):
                        def show_details():
                            details_dialog = QDialog(dialog)
                            details_dialog.setWindowTitle(f"Variable: {var_name}")
                            details_dialog.setMinimumWidth(500)
                            
                            details_layout = QVBoxLayout()
                            
                            # Variable name and type
                            name_label = QLabel(f"<b>{var_name}</b>")
                            details_layout.addWidget(name_label)
                            
                            type_label = QLabel(f"Type: <code>{var_obj.type}</code>")
                            type_label.setWordWrap(True)
                            details_layout.addWidget(type_label)
                            
                            # Description
                            desc_label = QLabel(f"<br><b>Description:</b><br>{var_obj.desc}")
                            desc_label.setWordWrap(True)
                            desc_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
                            details_layout.addWidget(desc_label)
                            
                            # Close button
                            close_btn = QDialogButtonBox(QDialogButtonBox.Ok)
                            close_btn.accepted.connect(details_dialog.accept)
                            details_layout.addWidget(close_btn)
                            
                            details_dialog.setLayout(details_layout)
                            details_dialog.exec_()
                        return show_details
                    
                    input_details_btn.clicked.connect(make_input_details_handler(inp, wsvar))
                    input_row_layout.addWidget(input_details_btn)
            
            input_row.setLayout(input_row_layout)
            input_layout.addRow(f"{inp}:", input_row)
        
        input_group.setLayout(input_layout)
        scroll_layout.addWidget(input_group)

    # --- Generic Inputs (editable) ---
    gin_values = {}  # Store current values for GIN parameters
    gin_selected_types = {}  # Store selected type for multi-type GINs
    gin_editors = {}  # Store editor widgets
    
    if method_obj.gin:
        gin_group = QGroupBox("Generic Inputs")
        gin_layout = QFormLayout()
        
        for i, gin_name in enumerate(method_obj.gin):
            gin_type_str = method_obj.gin_type[i]
            gin_default = method_obj.gin_value[i]
            gin_desc = method_obj.gin_desc[i] if i < len(method_obj.gin_desc) else ""
            
            # Parse gin_type: could be single type, comma-separated types, or "Any"
            if gin_type_str == "Any":
                # Any means all workspace groups
                available_types = sorted(cxx.globals.workspace_groups())
            elif "," in gin_type_str:
                # Multiple types specified
                available_types = [t.strip() for t in gin_type_str.split(",")]
            else:
                # Single type
                available_types = [gin_type_str]
            
            # Determine initial type and value
            if gin_default is None:
                # No default - use first available type
                selected_type = available_types[0]
                gin_selected_types[gin_name] = selected_type
                # Try to create a default instance of the type
                try:
                    if hasattr(cxx, selected_type):
                        gin_values[gin_name] = getattr(cxx, selected_type)()
                    else:
                        gin_values[gin_name] = None
                except Exception:
                    gin_values[gin_name] = None
            else:
                # Extract actual value from Wsv wrapper
                if hasattr(gin_default, 'value'):
                    actual_value = gin_default.value
                    gin_values[gin_name] = actual_value
                    # Determine which type was used
                    actual_type = type(actual_value).__name__
                    if actual_type in available_types:
                        gin_selected_types[gin_name] = actual_type
                    else:
                        gin_selected_types[gin_name] = available_types[0]
                else:
                    gin_values[gin_name] = gin_default
                    gin_selected_types[gin_name] = available_types[0]
            
            # Create row with type selector (if multiple types) and edit button
            row_widget = QWidget()
            row_layout = QHBoxLayout()
            row_layout.setContentsMargins(0, 0, 0, 0)
            
            def format_gin_value(val):
                """Format a GIN value for display."""
                if val is None:
                    return "<None>"
                type_name = type(val).__name__
                # Show a preview of the value if it's simple
                try:
                    val_str = str(val)
                    if len(val_str) > 50:
                        val_str = val_str[:50] + "..."
                    return f"{val_str}  ({type_name})"
                except Exception:
                    return f"<{type_name}>"
            
            # Type selector dropdown (if multiple types available)
            from PyQt5.QtWidgets import QComboBox
            type_combo = None
            if len(available_types) > 1:
                type_combo = QComboBox()
                type_combo.addItems(available_types)
                type_combo.setCurrentText(gin_selected_types[gin_name])
                type_combo.setMaximumWidth(200)
                row_layout.addWidget(type_combo)
                
                # Handler for type change - reinitialize value with new type
                def make_type_change_handler(gname, vlabel):
                    def on_type_changed(new_type):
                        gin_selected_types[gname] = new_type
                        # Create new instance of selected type
                        try:
                            if hasattr(cxx, new_type):
                                gin_values[gname] = getattr(cxx, new_type)()
                                vlabel.setText(format_gin_value(gin_values[gname]))
                            else:
                                gin_values[gname] = None
                                vlabel.setText("<Type not available>")
                        except Exception as e:
                            gin_values[gname] = None
                            vlabel.setText(f"<Error: {e}>")
                    return on_type_changed
            
            value_label = QLabel()
            value_label.setText(format_gin_value(gin_values[gin_name]))
            value_label.setStyleSheet("color: #555;")
            
            # Connect type change handler if we have a combo
            if type_combo is not None:
                type_combo.currentTextChanged.connect(make_type_change_handler(gin_name, value_label))
            
            edit_btn = QPushButton("Edit...")
            edit_btn.setMaximumWidth(80)
            
            def make_edit_callback(gname, avail_types, vlabel, tcombo):
                def edit_gin():
                    # Get current selected type
                    if tcombo is not None:
                        selected_type = tcombo.currentText()
                        gin_selected_types[gname] = selected_type
                    else:
                        selected_type = avail_types[0]
                    
                    current = gin_values[gname]
                    # If value is None or type changed, create new instance
                    if current is None or (tcombo and type(current).__name__ != selected_type):
                        try:
                            if hasattr(cxx, selected_type):
                                current = getattr(cxx, selected_type)()
                            else:
                                QMessageBox.warning(dialog, "Cannot edit",
                                                  f"Type '{selected_type}' not available")
                                return
                        except Exception as e:
                            QMessageBox.warning(dialog, "Cannot edit",
                                              f"Failed to initialize '{selected_type}': {e}")
                            return
                    
                    new_value = editors.edit(current, parent=dialog)
                    if new_value is not None:
                        gin_values[gname] = new_value
                        vlabel.setText(format_gin_value(new_value))
                return edit_gin
            
            edit_btn.clicked.connect(make_edit_callback(gin_name, available_types, value_label, type_combo))
            
            row_layout.addWidget(value_label, stretch=1)
            row_layout.addWidget(edit_btn)
            
            # Add "More details..." button if there's a description
            if gin_desc:
                details_btn = QPushButton("?")
                details_btn.setMaximumWidth(30)
                details_btn.setToolTip("Show parameter description")
                
                def make_details_handler(param_name, param_desc, param_types):
                    def show_details():
                        details_dialog = QDialog(dialog)
                        details_dialog.setWindowTitle(f"Parameter: {param_name}")
                        details_dialog.setMinimumWidth(500)
                        
                        details_layout = QVBoxLayout()
                        
                        # Parameter name and type
                        name_label = QLabel(f"<b>{param_name}</b>")
                        details_layout.addWidget(name_label)
                        
                        if len(param_types) == 1:
                            type_label = QLabel(f"Type: <code>{param_types[0]}</code>")
                        else:
                            type_label = QLabel(f"Types: <code>{', '.join(param_types)}</code>")
                        type_label.setWordWrap(True)
                        details_layout.addWidget(type_label)
                        
                        # Description
                        desc_label = QLabel(f"<br><b>Description:</b><br>{param_desc}")
                        desc_label.setWordWrap(True)
                        desc_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
                        details_layout.addWidget(desc_label)
                        
                        # Close button
                        close_btn = QDialogButtonBox(QDialogButtonBox.Ok)
                        close_btn.accepted.connect(details_dialog.accept)
                        details_layout.addWidget(close_btn)
                        
                        details_dialog.setLayout(details_layout)
                        details_dialog.exec_()
                    return show_details
                
                details_btn.clicked.connect(make_details_handler(gin_name, gin_desc, available_types))
                row_layout.addWidget(details_btn)
            
            row_widget.setLayout(row_layout)
            
            label_text = f"{gin_name}"
            if len(available_types) == 1:
                label_text += f" ({available_types[0]})"
            else:
                label_text += f" (multiple types)"
            
            gin_layout.addRow(label_text, row_widget)
            gin_editors[gin_name] = (value_label, edit_btn)
        
        gin_group.setLayout(gin_layout)
        scroll_layout.addWidget(gin_group)

    # --- Workspace Outputs (will be set) ---
    if method_obj.output:
        # Get workspace variable descriptions (reuse from inputs section)
        if 'wsvars' not in locals():
            wsvars = cxx.globals.workspace_variables()
        
        output_group = QGroupBox("Workspace Outputs (will be set)")
        output_layout = QFormLayout()
        for out in method_obj.output:
            # Row with status display and details button
            output_row = QWidget()
            output_row_layout = QHBoxLayout()
            output_row_layout.setContentsMargins(0, 0, 0, 0)
            
            # Check if it's also an input (in-place modification)
            if out in method_obj.input:
                label = QLabel("(modified in-place)")
            else:
                label = QLabel("(will be created/set)")
            label.setStyleSheet("color: #555;")
            output_row_layout.addWidget(label, stretch=1)
            
            # Add "?" button if variable has description
            if out in wsvars:
                wsvar = wsvars[out]
                if hasattr(wsvar, 'desc') and wsvar.desc:
                    output_details_btn = QPushButton("?")
                    output_details_btn.setMaximumWidth(30)
                    output_details_btn.setToolTip("Show variable description")
                    
                    def make_output_details_handler(var_name, var_obj):
                        def show_details():
                            details_dialog = QDialog(dialog)
                            details_dialog.setWindowTitle(f"Variable: {var_name}")
                            details_dialog.setMinimumWidth(500)
                            
                            details_layout = QVBoxLayout()
                            
                            # Variable name and type
                            name_label = QLabel(f"<b>{var_name}</b>")
                            details_layout.addWidget(name_label)
                            
                            type_label = QLabel(f"Type: <code>{var_obj.type}</code>")
                            type_label.setWordWrap(True)
                            details_layout.addWidget(type_label)
                            
                            # Description
                            desc_label = QLabel(f"<br><b>Description:</b><br>{var_obj.desc}")
                            desc_label.setWordWrap(True)
                            desc_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
                            details_layout.addWidget(desc_label)
                            
                            # Close button
                            close_btn = QDialogButtonBox(QDialogButtonBox.Ok)
                            close_btn.accepted.connect(details_dialog.accept)
                            details_layout.addWidget(close_btn)
                            
                            details_dialog.setLayout(details_layout)
                            details_dialog.exec_()
                        return show_details
                    
                    output_details_btn.clicked.connect(make_output_details_handler(out, wsvar))
                    output_row_layout.addWidget(output_details_btn)
            
            output_row.setLayout(output_row_layout)
            output_layout.addRow(f"{out}:", output_row)
        
        output_group.setLayout(output_layout)
        scroll_layout.addWidget(output_group)

    # --- Generic Outputs ---
    gout_names = {}  # Store names for GOUT parameters
    
    if method_obj.gout:
        gout_group = QGroupBox("Generic Outputs")
        gout_layout = QFormLayout()
        
        for i, gout_type in enumerate(method_obj.gout_type):
            gout_desc = method_obj.gout_desc[i] if i < len(method_obj.gout_desc) else ""
            default_name = method_obj.gout[i] if i < len(method_obj.gout) else f"gout_{i}"
            
            # Row with text field and details button
            gout_row = QWidget()
            gout_row_layout = QHBoxLayout()
            gout_row_layout.setContentsMargins(0, 0, 0, 0)
            
            # Text field for the output variable name
            name_edit = QLineEdit(default_name)
            gout_names[i] = name_edit
            gout_row_layout.addWidget(name_edit, stretch=1)
            
            # Add "More details..." button if there's a description
            if gout_desc:
                gout_details_btn = QPushButton("?")
                gout_details_btn.setMaximumWidth(30)
                gout_details_btn.setToolTip("Show output description")
                
                def make_gout_details_handler(param_name, param_type, param_desc):
                    def show_details():
                        details_dialog = QDialog(dialog)
                        details_dialog.setWindowTitle(f"Output: {param_name}")
                        details_dialog.setMinimumWidth(500)
                        
                        details_layout = QVBoxLayout()
                        
                        # Parameter name and type
                        name_label = QLabel(f"<b>{param_name}</b>")
                        details_layout.addWidget(name_label)
                        
                        type_label = QLabel(f"Type: <code>{param_type}</code>")
                        type_label.setWordWrap(True)
                        details_layout.addWidget(type_label)
                        
                        # Description
                        desc_label = QLabel(f"<br><b>Description:</b><br>{param_desc}")
                        desc_label.setWordWrap(True)
                        desc_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
                        details_layout.addWidget(desc_label)
                        
                        # Close button
                        close_btn = QDialogButtonBox(QDialogButtonBox.Ok)
                        close_btn.accepted.connect(details_dialog.accept)
                        details_layout.addWidget(close_btn)
                        
                        details_dialog.setLayout(details_layout)
                        details_dialog.exec_()
                    return show_details
                
                gout_details_btn.clicked.connect(make_gout_details_handler(default_name, gout_type, gout_desc))
                gout_row_layout.addWidget(gout_details_btn)
            
            gout_row.setLayout(gout_row_layout)
            
            label_text = f"{gout_type}"
            gout_layout.addRow(label_text, gout_row)
        
        gout_group.setLayout(gout_layout)
        scroll_layout.addWidget(gout_group)

    # Buttons
    button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    button_box.accepted.connect(lambda: execute_method())
    button_box.rejected.connect(dialog.reject)
    main.addWidget(button_box)

    dialog.setLayout(main)

    def execute_method():
        """Execute the method with current parameters."""
        try:
            # Get the method callable from workspace
            method_callable = getattr(ws, method_name)
            
            # Prepare keyword arguments for GIN
            kwargs = {}
            for i, gin_name in enumerate(method_obj.gin):
                if gin_name in gin_values:
                    value = gin_values[gin_name]
                    # Only pass non-None values, or if the method requires it
                    if value is not None or method_obj.gin_value[i] is None:
                        # Use the actual value if set, otherwise None will be passed
                        # (which may cause an error if it's required)
                        kwargs[gin_name] = value
            
            # Handle GOUT - these are output variable names
            # GOUT names are passed as keyword arguments with the GOUT parameter name
            if method_obj.gout:
                for i, gout_name in enumerate(method_obj.gout):
                    var_name = gout_names[i].text().strip()
                    if not var_name:
                        QMessageBox.warning(dialog, "Invalid GOUT",
                                          f"Generic output '{gout_name}' needs a variable name")
                        return
                    # GOUT is typically passed as keyword with gout parameter name = ws variable name
                    kwargs[gout_name] = var_name
            
            # Call the method
            method_callable(**kwargs)
            
            # Show success with info about outputs
            outputs_info = ""
            if method_obj.output:
                outputs_info = "\n\nModified workspace variables:\n" + "\n".join(f"  - {out}" for out in method_obj.output)
            if method_obj.gout:
                gout_vars = [gout_names[i].text().strip() for i in range(len(method_obj.gout))]
                outputs_info += "\n\nGeneric outputs:\n" + "\n".join(f"  - {v}" for v in gout_vars)
            
            QMessageBox.information(dialog, "Success",
                                  f"Method '{method_name}' executed successfully{outputs_info}")
            dialog.accept()
            
        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            
            # Show error with scrollable details
            msg_box = QMessageBox(dialog)
            msg_box.setIcon(QMessageBox.Critical)
            msg_box.setWindowTitle("Execution Failed")
            error_lines = str(e).split('\n')
            last_line = error_lines[-1].strip() if error_lines else str(e)
            msg_box.setText(f"Failed to execute '{method_name}':\n\n{last_line}")
            msg_box.setDetailedText(tb)
            msg_box.setStyleSheet("QTextEdit { min-width: 600px; min-height: 300px; }")
            msg_box.exec_()

    dialog.exec_()
