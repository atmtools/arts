"""Editor for ARTS option group enum values.

Provides a dropdown menu to select from available enum options.
"""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QDialogButtonBox, QApplication
)

__all__ = ["edit"]


def edit(value, parent=None):
    """
    Edit an ARTS option group enum value using a dropdown selector.

    Parameters
    ----------
    value : option enum
        An instance of an ARTS option group enum (e.g., ray_path_observer_agendaSetGeometricMaxStep)
    parent : QWidget, optional
        Parent widget for the dialog

    Returns
    -------
    option enum or None
        The selected option value if accepted, None if cancelled
    """
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    
    # Get the type and available options
    option_type = type(value)
    type_name = option_type.__name__
    
    try:
        options = option_type.get_options()
    except Exception as e:
        # Fallback: if get_options() doesn't work, can't edit
        from .common import create_description_dialog
        dialog = create_description_dialog(
            title="Cannot Edit Option",
            name=type_name,
            type_str=type_name,
            description=f"Unable to retrieve options: {e}"
        )
        dialog.exec_()
        return None
    
    # Create dialog
    dialog = QDialog(parent)
    dialog.setWindowTitle(f"Edit {type_name}")
    dialog.resize(500, 150)
    
    layout = QVBoxLayout()
    
    # Type label
    type_label = QLabel(f"<b>Type:</b> {type_name}")
    layout.addWidget(type_label)
    
    # Current value label
    current_label = QLabel(f"<b>Current value:</b> {str(value)}")
    layout.addWidget(current_label)
    
    # Dropdown for selection
    selector_layout = QHBoxLayout()
    selector_label = QLabel("Select option:")
    selector_layout.addWidget(selector_label)
    
    combo = QComboBox()
    combo.setMaxVisibleItems(20)  # Make it scrollable for long lists
    
    # Populate combo box with options
    option_strings = []
    current_index = 0
    for i, opt in enumerate(options):
        opt_str = str(opt)
        option_strings.append(opt_str)
        combo.addItem(opt_str)
        
        # Find current selection
        if str(opt) == str(value):
            current_index = i
    
    combo.setCurrentIndex(current_index)
    selector_layout.addWidget(combo, stretch=1)
    layout.addLayout(selector_layout)
    
    # Add some spacing
    layout.addStretch()
    
    # Buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    # Show dialog and get result
    if dialog.exec_() == QDialog.Accepted:
        selected_index = combo.currentIndex()
        return options[selected_index]
    
    return None
