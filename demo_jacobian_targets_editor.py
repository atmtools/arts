#!/usr/bin/env python3
"""Interactive demo for the JacobianTargets editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class JacobianTargetsEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("JacobianTargets Editor Demo")
        self.setMinimumSize(700, 500)
        
        # Create JacobianTargets instance
        self.targets = arts.JacobianTargets()
        # All arrays start empty - user can add targets
        
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
        btn = QPushButton("Edit JacobianTargets")
        btn.clicked.connect(self.edit_targets)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current values."""
        text = "Current JacobianTargets:\n"
        text += "=" * 50 + "\n\n"
        text += "Jacobian targets for retrieval/sensitivity analysis:\n\n"
        text += f"Atmospheric targets:    {len(self.targets.atm)} targets\n"
        text += f"Surface targets:        {len(self.targets.surf)} targets\n"
        text += f"Subsurface targets:     {len(self.targets.subsurf)} targets\n"
        text += f"Line targets:           {len(self.targets.line)} targets\n"
        text += f"Sensor targets:         {len(self.targets.sensor)} targets\n"
        text += f"Error targets:          {len(self.targets.error)} targets\n"
        text += "\n" + "=" * 50 + "\n"
        text += "\nClick 'Edit JacobianTargets' to open the editor.\n"
        text += "Double-click any target array to add/edit targets.\n"
        text += "\nExample targets:\n"
        text += "  - Atmospheric: temperature, gas concentrations\n"
        text += "  - Surface: emissivity, reflectivity\n"
        text += "  - Line: line strength, broadening"
        self.text.setText(text)
    
    def edit_targets(self):
        """Open the JacobianTargets editor."""
        result = edit(self.targets, parent=self)
        if result is not None:
            self.targets = result
            self.update_display()
            print("\n✓ JacobianTargets updated!")
            total = (len(self.targets.atm) + len(self.targets.surf) + 
                    len(self.targets.subsurf) + len(self.targets.line) +
                    len(self.targets.sensor) + len(self.targets.error))
            print(f"  Total targets: {total}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = JacobianTargetsEditorDemo()
    demo.show()
    sys.exit(app.exec_())
