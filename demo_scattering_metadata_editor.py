#!/usr/bin/env python3
"""Interactive demo for the ScatteringMetaData editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class ScatteringMetaDataEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ScatteringMetaData Editor Demo")
        self.setMinimumSize(700, 500)
        
        # Create ScatteringMetaData with sample data
        self.metadata = arts.ScatteringMetaData()
        self.metadata.description = "Ice particle"
        self.metadata.source = "Laboratory measurements"
        self.metadata.refr_index = "1.33 + 0.001i"
        self.metadata.mass = 1.5e-10  # kg
        self.metadata.diameter_max = 500e-6  # meters
        self.metadata.diameter_volume_equ = 450e-6  # meters
        self.metadata.diameter_area_equ_aerodynamical = 480e-6  # meters
        
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
        btn = QPushButton("Edit ScatteringMetaData")
        btn.clicked.connect(self.edit_metadata)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current values."""
        text = "Current ScatteringMetaData:\n"
        text += "=" * 60 + "\n\n"
        text += f"Description:    {self.metadata.description}\n"
        text += f"Source:         {self.metadata.source}\n"
        text += f"Refr. Index:    {self.metadata.refr_index}\n"
        text += f"Mass:           {self.metadata.mass:.3e} kg\n"
        text += f"Max Diameter:   {self.metadata.diameter_max*1e6:.1f} µm\n"
        text += f"Volume Equiv.:  {self.metadata.diameter_volume_equ*1e6:.1f} µm\n"
        text += f"Area Equiv.:    {self.metadata.diameter_area_equ_aerodynamical*1e6:.1f} µm\n"
        text += "\n" + "=" * 60 + "\n"
        text += "\nMetadata for scattering particles used in radiative transfer.\n"
        text += "Click 'Edit ScatteringMetaData' to modify properties.\n"
        text += "Double-click any Value to edit."
        self.text.setText(text)
    
    def edit_metadata(self):
        """Open the ScatteringMetaData editor."""
        result = edit(self.metadata, parent=self)
        if result is not None:
            self.metadata = result
            self.update_display()
            print("\n✓ ScatteringMetaData updated!")
            print(f"  Description: {self.metadata.description}")
            print(f"  Mass: {self.metadata.mass:.3e} kg")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = ScatteringMetaDataEditorDemo()
    demo.show()
    sys.exit(app.exec_())
