#!/usr/bin/env python3
"""Interactive demo for the SpeciesTag editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class SpeciesTagEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SpeciesTag Editor Demo")
        self.setMinimumSize(600, 400)
        
        # Create SpeciesTag with sample data
        # SpeciesTag properties are read-only, create from string
        self.tag = arts.SpeciesTag("H2O")
        
        # Create UI
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        
        # Text display
        self.text = QTextEdit()
        self.text.setReadOnly(True)
        self.update_display()
        layout.addWidget(self.text)
        
        # Edit button
        btn = QPushButton("Edit SpeciesTag")
        btn.clicked.connect(self.edit_tag)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current values."""
        text = "Current SpeciesTag:\n"
        text += "=" * 50 + "\n\n"
        text += f"Species:          {self.tag.spec}\n"
        text += f"Type:             {self.tag.type}\n"
        text += f"CIA 2nd Species:  {self.tag.cia_2nd_species}\n"
        text += f"Species Index:    {self.tag.spec_ind}\n"
        text += f"Full Name:        {self.tag.full_name}\n"
        text += "\n" + "=" * 50 + "\n"
        text += "\nClick 'Edit SpeciesTag' to open the editor.\n"
        text += "Try changing the species (H2O, O2, N2, etc.)\n"
        text += "or the tag type (Plain, Predefined, CIA, etc.)."
        self.text.setText(text)
    
    def edit_tag(self):
        """Open the SpeciesTag editor."""
        result = edit(self.tag, parent=self)
        if result is not None:
            self.tag = result
            self.update_display()
            print("\n✓ SpeciesTag updated!")
            print(f"  Species: {self.tag.spec}")
            print(f"  Full name: {self.tag.full_name}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = SpeciesTagEditorDemo()
    demo.show()
    sys.exit(app.exec_())
