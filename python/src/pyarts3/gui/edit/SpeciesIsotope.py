"""Editor for SpeciesIsotope type.

SpeciesIsotope can be selected from all_isotopologues().
It has read-only fields:
- spec: Species
- isotname: String
- mass: Numeric
- gi: Index
"""

from PyQt5.QtWidgets import (
    QApplication, QDialog, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QDialogButtonBox
)

__all__ = ["edit"]


def edit(value, parent=None):
    """
    Edit a SpeciesIsotope value.
    
    Parameters
    ----------
    value : SpeciesIsotope
        The species isotope to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    SpeciesIsotope or None
        The edited value if accepted, None if cancelled
    """
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit SpeciesIsotope")
    dialog.resize(500, 150)
    
    layout = QVBoxLayout()
    
    # Info
    info = QLabel(f"Current: {value.spec} isotope {value.isotname} (mass: {value.mass}, gi: {value.gi})")
    info.setStyleSheet("color: #555; font-weight: bold;")
    layout.addWidget(info)
    
    # Get all available isotopologues
    from pyarts3 import arts
    all_isotopologues = arts.globals.all_isotopologues()
    
    # Dropdown for selection
    combo_layout = QHBoxLayout()
    combo_layout.addWidget(QLabel("Select Isotopologue:"))
    combo = QComboBox()
    combo.setMaxVisibleItems(20)  # Make it scrollable
    
    # Populate dropdown
    current_index = 0
    for i, isot in enumerate(all_isotopologues):
        combo.addItem(str(isot), isot)  # Display text and store the object
        if str(isot) == str(value):
            current_index = i
    
    combo.setCurrentIndex(current_index)
    combo_layout.addWidget(combo)
    layout.addLayout(combo_layout)
    
    # Help text
    help_text = QLabel(f"Choose from {len(all_isotopologues)} available isotopologues")
    help_text.setStyleSheet("color: #888; font-size: 10px;")
    layout.addWidget(help_text)
    
    layout.addStretch()
    
    # Dialog buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    if dialog.exec_() == QDialog.Accepted:
        return combo.currentData()  # Return the stored SpeciesIsotope object
    return None
