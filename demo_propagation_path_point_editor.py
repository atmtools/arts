#!/usr/bin/env python3
"""Interactive demo for the PropagationPathPoint editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class PropagationPathPointEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PropagationPathPoint Editor Demo")
        self.setMinimumSize(650, 450)
        
        # Create PropagationPathPoint with sample data
        self.ppp = arts.PropagationPathPoint()
        self.ppp.pos = arts.Vector3([6.371e6, 0.0, 0.0])  # Earth radius at equator
        self.ppp.pos_type = arts.PathPositionType.atm
        self.ppp.los = arts.Vector2([90.0, 0.0])  # Zenith, azimuth
        self.ppp.los_type = arts.PathPositionType.atm
        self.ppp.nreal = 1.0003  # Refractive index
        self.ppp.ngroup = 1.0003
        
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
        btn = QPushButton("Edit PropagationPathPoint")
        btn.clicked.connect(self.edit_ppp)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current values."""
        text = "Current PropagationPathPoint:\n"
        text += "=" * 50 + "\n\n"
        text += f"Position:       {self.ppp.pos}\n"
        text += f"Position Type:  {self.ppp.pos_type}\n"
        text += f"Line-of-sight:  {self.ppp.los}\n"
        text += f"LOS Type:       {self.ppp.los_type}\n"
        text += f"n_real:         {self.ppp.nreal:.6f}\n"
        text += f"n_group:        {self.ppp.ngroup:.6f}\n"
        text += "\n" + "=" * 50 + "\n"
        text += "\nClick 'Edit PropagationPathPoint' to open the editor.\n"
        text += "Double-click any Value to edit (Vector3, Vector2, enums, etc.)."
        self.text.setText(text)
    
    def edit_ppp(self):
        """Open the PropagationPathPoint editor."""
        result = edit(self.ppp, parent=self)
        if result is not None:
            self.ppp = result
            self.update_display()
            print("\n✓ PropagationPathPoint updated!")
            print(f"  Position: {self.ppp.pos}")
            print(f"  LOS: {self.ppp.los}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = PropagationPathPointEditorDemo()
    demo.show()
    sys.exit(app.exec_())
