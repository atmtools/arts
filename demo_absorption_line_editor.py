#!/usr/bin/env python3
"""Interactive demo for the AbsorptionLine editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class AbsorptionLineEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("AbsorptionLine Editor Demo")
        self.setMinimumSize(700, 550)
        
        # Create AbsorptionLine with sample spectroscopic data
        self.line = arts.AbsorptionLine()
        self.line.f0 = 183.31012e9  # H2O line at 183.31 GHz
        self.line.a = 1.4e-20  # Line strength
        self.line.e0 = 136.76  # Lower state energy (K)
        self.line.gu = 16.0  # Upper state degeneracy
        self.line.gl = 14.0  # Lower state degeneracy
        # ls (LineShapeModel), z (ZeemanLineModel), qn (QuantumState) have defaults
        
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
        btn = QPushButton("Edit AbsorptionLine")
        btn.clicked.connect(self.edit_line)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current values."""
        text = "Current AbsorptionLine:\n"
        text += "=" * 50 + "\n\n"
        text += "Basic Spectroscopic Parameters:\n"
        text += f"  f0 (center freq):     {self.line.f0:.6e} Hz\n"
        text += f"                        ({self.line.f0/1e9:.6f} GHz)\n"
        text += f"  a (line strength):    {self.line.a:.6e}\n"
        text += f"  e0 (lower energy):    {self.line.e0:.2f} K\n"
        text += f"  gu (upper degenr):    {self.line.gu:.1f}\n"
        text += f"  gl (lower degenr):    {self.line.gl:.1f}\n\n"
        text += "Advanced Parameters:\n"
        text += f"  ls (line shape):      {type(self.line.ls).__name__}\n"
        text += f"  z (Zeeman model):     {type(self.line.z).__name__}\n"
        text += f"  qn (quantum state):   {type(self.line.qn).__name__}\n"
        text += "\n" + "=" * 50 + "\n"
        text += "\nClick 'Edit AbsorptionLine' to open the editor.\n"
        text += "This editor allows you to modify all spectroscopic\n"
        text += "parameters of an individual absorption line.\n\n"
        text += "Double-click ls/z/qn to edit complex line shape,\n"
        text += "Zeeman splitting, or quantum number assignments."
        self.text.setText(text)
    
    def edit_line(self):
        """Open the AbsorptionLine editor."""
        result = edit(self.line, parent=self)
        if result is not None:
            self.line = result
            self.update_display()
            print("\n✓ AbsorptionLine updated!")
            print(f"  Frequency: {self.line.f0/1e9:.6f} GHz")
            print(f"  Strength: {self.line.a:.6e}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = AbsorptionLineEditorDemo()
    demo.show()
    sys.exit(app.exec_())
