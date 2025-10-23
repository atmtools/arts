#!/usr/bin/env python3
"""Interactive demo for the XsecRecord editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class XsecRecordEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("XsecRecord Editor Demo")
        self.setMinimumSize(750, 550)
        
        # Create XsecRecord with sample data
        # XsecRecord properties are read-only after creation
        self.xsec = arts.XsecRecord()
        # Properties are initialized to defaults and can only be viewed
        
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
        btn = QPushButton("Edit XsecRecord")
        btn.clicked.connect(self.edit_xsec)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current values."""
        text = "Current XsecRecord:\n"
        text += "=" * 60 + "\n\n"
        text += f"Species:                {self.xsec.species}\n"
        text += f"Version:                {self.xsec.version}\n"
        text += f"Min Pressures:          {len(self.xsec.fitminpressures)} values\n"
        text += f"Max Pressures:          {len(self.xsec.fitmaxpressures)} values\n"
        text += f"Min Temperatures:       {len(self.xsec.fitmintemperatures)} values\n"
        text += f"Max Temperatures:       {len(self.xsec.fitmaxtemperatures)} values\n"
        text += f"Fit Coefficients:       {len(self.xsec.fitcoeffs)} fields\n"
        text += "\n" + "=" * 60 + "\n"
        text += "\nCross-Section Absorption Record\n\n"
        text += "Cross-section data provides pre-calculated absorption\n"
        text += "coefficients with polynomial fit coefficients as a function\n"
        text += "of temperature and pressure.\n\n"
        text += "This is used for species with complex absorption spectra\n"
        text += "where line-by-line calculations would be too expensive.\n\n"
        text += "Common species:\n"
        text += "  • CFC-11, CFC-12: Chlorofluorocarbons\n"
        text += "  • O3: Ozone (UV-Visible bands)\n"
        text += "  • N2O5, NO2: Nitrogen oxides\n\n"
        text += "NOTE: XsecRecord properties are read-only after creation.\n"
        text += "Click 'Edit XsecRecord' to view the structure.\n"
        text += "In practice, XsecRecord objects are loaded from files."
        self.text.setText(text)
    
    def edit_xsec(self):
        """Open the XsecRecord editor."""
        result = edit(self.xsec, parent=self)
        if result is not None:
            self.xsec = result
            self.update_display()
            print("\n✓ XsecRecord updated!")
            print(f"  Species: {self.xsec.species}")
            print(f"  Fit coefficients: {len(self.xsec.fitcoeffs)}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = XsecRecordEditorDemo()
    demo.show()
    sys.exit(app.exec_())
