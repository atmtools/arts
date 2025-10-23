#!/usr/bin/env python3
"""Interactive demo for the Sun editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class SunEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Sun Editor Demo")
        self.setMinimumSize(600, 400)
        
        # Create Sun instance with some sample data
        self.sun = arts.Sun()
        self.sun.description = "Sun for radiative transfer"
        self.sun.distance = 1.496e11  # AU in meters
        self.sun.latitude = 0.0
        self.sun.longitude = 0.0
        self.sun.radius = 6.96e8  # meters
        # spectrum is a Matrix - start empty
        
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
        btn = QPushButton("Edit Sun")
        btn.clicked.connect(self.edit_sun)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current Sun values."""
        text = "Current Sun Properties:\n"
        text += "=" * 50 + "\n\n"
        text += f"Description:  {self.sun.description}\n"
        text += f"Distance:     {self.sun.distance:.3e} m\n"
        text += f"Latitude:     {self.sun.latitude:.2f}°\n"
        text += f"Longitude:    {self.sun.longitude:.2f}°\n"
        text += f"Radius:       {self.sun.radius:.3e} m\n"
        text += f"Spectrum:     {self.sun.spectrum.nrows}x{self.sun.spectrum.ncols} Matrix\n"
        text += "\n" + "=" * 50 + "\n"
        text += "\nClick 'Edit Sun' to open the editor dialog.\n"
        text += "Double-click any Value to edit nested objects."
        self.text.setText(text)
    
    def edit_sun(self):
        """Open the Sun editor."""
        result = edit(self.sun, parent=self)
        if result is not None:
            self.sun = result
            self.update_display()
            print("\n✓ Sun updated!")
            print(f"  Distance: {self.sun.distance:.3e} m")
            print(f"  Position: ({self.sun.latitude:.2f}°, {self.sun.longitude:.2f}°)")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = SunEditorDemo()
    demo.show()
    sys.exit(app.exec_())
