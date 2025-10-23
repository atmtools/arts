#!/usr/bin/env python3
"""Interactive demo for the SensorObsel editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class SensorObselEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SensorObsel Editor Demo")
        self.setMinimumSize(650, 400)
        
        # Create SensorObsel with sample data
        self.obsel = arts.SensorObsel()
        # Note: SensorObsel fields are read-only after creation
        # They start with default empty values
        # User can edit them through the editor dialog
        
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
        btn = QPushButton("Edit SensorObsel")
        btn.clicked.connect(self.edit_obsel)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current values."""
        text = "Current SensorObsel:\n"
        text += "=" * 50 + "\n\n"
        text += f"Frequency Grid:  {len(self.obsel.f_grid)} points\n"
        if len(self.obsel.f_grid) > 0:
            text += f"  Range: {self.obsel.f_grid[0]:.3e} - {self.obsel.f_grid[-1]:.3e} Hz\n"
        text += f"PosLos Vector:   {len(self.obsel.poslos)} elements\n"
        text += f"Weight Matrix:   {self.obsel.weight_matrix.shape}\n"
        text += "\n" + "=" * 50 + "\n"
        text += "\nClick 'Edit SensorObsel' to open the editor.\n"
        text += "Double-click f_grid to edit the frequency grid,\n"
        text += "or poslos/weight_matrix to configure sensor geometry."
        self.text.setText(text)
    
    def edit_obsel(self):
        """Open the SensorObsel editor."""
        result = edit(self.obsel, parent=self)
        if result is not None:
            self.obsel = result
            self.update_display()
            print("\n✓ SensorObsel updated!")
            print(f"  Frequency points: {len(self.obsel.f_grid)}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = SensorObselEditorDemo()
    demo.show()
    sys.exit(app.exec_())
