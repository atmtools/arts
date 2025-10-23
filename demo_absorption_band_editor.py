#!/usr/bin/env python3
"""Interactive demo for the AbsorptionBand editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class AbsorptionBandEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("AbsorptionBand Editor Demo")
        self.setMinimumSize(650, 450)
        
        # Create AbsorptionBand with sample configuration
        self.band = arts.AbsorptionBand()
        self.band.cutoff = arts.LineByLineCutoffType.None_
        self.band.cutoff_value = float('inf')
        # lines starts as empty array
        self.band.lineshape = arts.LineByLineLineshape()
        
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
        btn = QPushButton("Edit AbsorptionBand")
        btn.clicked.connect(self.edit_band)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current values."""
        text = "Current AbsorptionBand:\n"
        text += "=" * 50 + "\n\n"
        text += f"Cutoff Type:    {self.band.cutoff}\n"
        text += f"Cutoff Value:   {self.band.cutoff_value}\n"
        text += f"Lines:          {len(self.band.lines)} absorption lines\n"
        text += f"Line Shape:     {self.band.lineshape}\n"
        text += "\n" + "=" * 50 + "\n"
        text += "\nClick 'Edit AbsorptionBand' to open the editor.\n"
        text += "Double-click 'lines' to add/edit individual absorption lines.\n"
        text += "Each line contains spectroscopic parameters (f0, a, e0, etc.)\n\n"
        text += "Cutoff types control line wing calculations:\n"
        text += "  - None: No cutoff\n"
        text += "  - ByLine: Per-line cutoff\n"
        text += "  - ByBand: Band-wide cutoff"
        self.text.setText(text)
    
    def edit_band(self):
        """Open the AbsorptionBand editor."""
        result = edit(self.band, parent=self)
        if result is not None:
            self.band = result
            self.update_display()
            print("\n✓ AbsorptionBand updated!")
            print(f"  Number of lines: {len(self.band.lines)}")
            print(f"  Cutoff: {self.band.cutoff}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = AbsorptionBandEditorDemo()
    demo.show()
    sys.exit(app.exec_())
