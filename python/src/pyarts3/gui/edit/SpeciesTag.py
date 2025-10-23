"""Editor for SpeciesTag type."""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QLineEdit,
    QDialogButtonBox, QFormLayout, QTextEdit
)


def edit(value, parent=None):
    """
    Edit a SpeciesTag object.
    
    SpeciesTag is constructed from a string specification that determines
    all its properties. The editor allows editing the full_name string,
    which is then parsed to create a new SpeciesTag.
    
    Examples:
    - "H2O" - Plain water vapor
    - "O2-66" - Oxygen isotopologue
    - "H2O-CIA-H2O" - Collision-induced absorption
    - "O3-*-183000000000" - Predefined frequency
    
    Parameters
    ----------
    value : SpeciesTag
        The SpeciesTag object to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    SpeciesTag or None
        Modified SpeciesTag if OK clicked, None if cancelled
    """
    from pyarts3 import arts
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit SpeciesTag")
    dialog.setMinimumWidth(600)
    dialog.setMinimumHeight(400)
    
    layout = QVBoxLayout(dialog)
    
    # Current values display
    info = QTextEdit()
    info.setReadOnly(True)
    info.setMaximumHeight(150)
    current_info = f"""Current SpeciesTag:
  Species:          {value.spec}
  Type:             {value.type}
  CIA 2nd Species:  {value.cia_2nd_species}
  Species Index:    {value.spec_ind}
  Full Name:        {value.full_name}

Note: SpeciesTag is constructed from the full_name string below."""
    info.setText(current_info)
    layout.addWidget(info)
    
    # Edit form
    form = QFormLayout()
    
    full_name_edit = QLineEdit(value.full_name)
    full_name_edit.setToolTip("Enter species tag string (e.g., 'H2O', 'O2-66', 'H2O-CIA-N2')")
    form.addRow("Full Name:", full_name_edit)
    
    layout.addLayout(form)
    
    # Help text
    help_text = QLabel(
        "Examples:\n"
        "  • H2O - Water vapor\n"
        "  • O2-66 - Oxygen isotopologue\n"
        "  • H2O-CIA-H2O - Collision-induced absorption\n"
        "  • O3-*-183000000000 - Predefined frequency"
    )
    help_text.setStyleSheet("QLabel { color: gray; font-size: 10pt; }")
    layout.addWidget(help_text)
    
    # OK/Cancel buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    if dialog.exec_() == QDialog.Accepted:
        try:
            # SpeciesTag is constructed from the full_name string
            result = arts.SpeciesTag(full_name_edit.text())
            return result
        except Exception as e:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.warning(
                dialog,
                "Invalid SpeciesTag",
                f"Cannot create SpeciesTag from '{full_name_edit.text()}':\n\n{str(e)}"
            )
            return None
    
    return None
