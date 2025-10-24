"""Editor for ARTS SpeciesIsotope values.

Provides a dropdown menu to select from available predefined species isotopes.
"""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QDialogButtonBox, QApplication, QLineEdit
)

__all__ = ["edit"]


def edit(value, parent=None):
    """
    Edit an ARTS SpeciesIsotope value using a dropdown selector with filter.

    Parameters
    ----------
    value : SpeciesIsotope
        An instance of SpeciesIsotope
    parent : QWidget, optional
        Parent widget for the dialog

    Returns
    -------
    SpeciesIsotope or None
        The selected SpeciesIsotope value if accepted, None if cancelled
    """
    import pyarts3.arts as arts
    
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    
    # Get predefined species isotopes
    predefined = arts.globals.all_isotopologues()
    
    # Create dialog
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit SpeciesIsotope")
    dialog.resize(500, 200)
    
    layout = QVBoxLayout()
    
    # Type label
    type_label = QLabel("<b>Type:</b> SpeciesIsotope")
    layout.addWidget(type_label)
    
    # Current value label
    current_label = QLabel(f"<b>Current value:</b> {str(value)}")
    layout.addWidget(current_label)
    
    # Filter box
    filter_layout = QHBoxLayout()
    filter_label = QLabel("Filter:")
    filter_layout.addWidget(filter_label)
    filter_box = QLineEdit()
    filter_box.setPlaceholderText("Type to filter isotopes...")
    filter_layout.addWidget(filter_box)
    layout.addLayout(filter_layout)
    
    # Dropdown for selection
    selector_layout = QHBoxLayout()
    selector_label = QLabel("Select isotope:")
    selector_layout.addWidget(selector_label)
    
    combo = QComboBox()
    combo.setMaxVisibleItems(20)  # Make it scrollable for long lists
    
    # Populate combo box with predefined isotopes
    all_isotopes = []
    current_index = 0
    for i, iso in enumerate(predefined):
        iso_str = str(iso)
        all_isotopes.append((iso_str, iso))
        if iso == value:
            current_index = i
    
    def update_combo(filter_text=""):
        combo.clear()
        new_index = 0
        filter_lower = filter_text.lower()
        for i, (iso_str, iso) in enumerate(all_isotopes):
            if not filter_text or filter_lower in iso_str.lower():
                combo.addItem(iso_str, userData=iso)
                if iso == value:
                    new_index = combo.count() - 1
        combo.setCurrentIndex(new_index)
    
    # Initial population
    update_combo()
    
    # Connect filter
    filter_box.textChanged.connect(update_combo)
    
    selector_layout.addWidget(combo)
    layout.addLayout(selector_layout)
    
    # Dialog buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    # Execute dialog
    if dialog.exec_() == QDialog.Accepted:
        selected_iso = combo.currentData()
        if selected_iso is not None:
            return selected_iso
    
    return None
