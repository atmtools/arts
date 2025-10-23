#!/usr/bin/env python3
"""Interactive demo for the ZeemanLineModel editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class ZeemanLineModelEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ZeemanLineModel Editor Demo")
        self.setMinimumSize(600, 400)
        
        # Create ZeemanLineModel with sample data
        self.zeeman = arts.ZeemanLineModel()
        self.zeeman.on = True
        self.zeeman.gu = 2.002  # Upper state Landé g-factor
        self.zeeman.gl = 1.998  # Lower state Landé g-factor
        
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
        btn = QPushButton("Edit ZeemanLineModel")
        btn.clicked.connect(self.edit_zeeman)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current values."""
        text = "Current ZeemanLineModel:\n"
        text += "=" * 60 + "\n\n"
        text += f"Zeeman Effect:  {'ON' if self.zeeman.on else 'OFF'}\n"
        text += f"Upper g-factor: {self.zeeman.gu:.6f}\n"
        text += f"Lower g-factor: {self.zeeman.gl:.6f}\n"
        text += "\n" + "=" * 60 + "\n"
        text += "\nThe Zeeman effect causes spectral line splitting in\n"
        text += "the presence of a magnetic field.\n\n"
        text += "The Landé g-factors determine the magnitude of the splitting:\n"
        text += "  • gu: Upper state g-factor\n"
        text += "  • gl: Lower state g-factor\n\n"
        text += "Click 'Edit ZeemanLineModel' to modify parameters.\n"
        text += "Double-click any Value to edit."
        self.text.setText(text)
    
    def edit_zeeman(self):
        """Open the ZeemanLineModel editor."""
        result = edit(self.zeeman, parent=self)
        if result is not None:
            self.zeeman = result
            self.update_display()
            print("\n✓ ZeemanLineModel updated!")
            print(f"  Zeeman effect: {'ON' if self.zeeman.on else 'OFF'}")
            print(f"  gu = {self.zeeman.gu:.6f}, gl = {self.zeeman.gl:.6f}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = ZeemanLineModelEditorDemo()
    demo.show()
    sys.exit(app.exec_())
