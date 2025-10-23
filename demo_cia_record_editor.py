#!/usr/bin/env python3
"""Interactive demo for the CIARecord editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class CIARecordEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CIARecord Editor Demo")
        self.setMinimumSize(700, 450)
        
        # Create CIARecord with sample data
        # CIARecord(data, species1, species2)
        data = arts.ArrayOfGriddedField2()  # Empty data array
        self.cia = arts.CIARecord(data, arts.SpeciesEnum.H2O, arts.SpeciesEnum.N2)
        
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
        btn = QPushButton("Edit CIARecord")
        btn.clicked.connect(self.edit_cia)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current values."""
        text = "Current CIARecord:\n"
        text += "=" * 60 + "\n\n"
        text += f"Species Pair:  {str(self.cia.specs[0])} + {str(self.cia.specs[1])}\n"
        text += f"Data Fields:   {len(self.cia.data)} GriddedField2 objects\n"
        text += "\n" + "=" * 60 + "\n"
        text += "\nCollision-Induced Absorption (CIA) Record\n\n"
        text += "CIA occurs when two molecules collide and temporarily form\n"
        text += "a complex that can absorb radiation at frequencies where\n"
        text += "neither molecule alone would absorb.\n\n"
        text += "Common pairs:\n"
        text += "  • N2-N2: Nitrogen dimer\n"
        text += "  • H2O-N2: Water-nitrogen complex\n"
        text += "  • O2-O2: Oxygen dimer\n"
        text += "  • CO2-CO2: Carbon dioxide dimer\n\n"
        text += "Click 'Edit CIARecord' to modify the record.\n"
        text += "Double-click 'specs' to edit the species pair.\n"
        text += "Double-click 'data' to add/edit absorption data."
        self.text.setText(text)
    
    def edit_cia(self):
        """Open the CIARecord editor."""
        result = edit(self.cia, parent=self)
        if result is not None:
            self.cia = result
            self.update_display()
            print("\n✓ CIARecord updated!")
            print(f"  Species: {str(self.cia.specs[0])} + {str(self.cia.specs[1])}")
            print(f"  Data fields: {len(self.cia.data)}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = CIARecordEditorDemo()
    demo.show()
    sys.exit(app.exec_())
