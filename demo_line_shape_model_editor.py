#!/usr/bin/env python3
"""Interactive demo for the LineShapeModel editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class LineShapeModelEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LineShapeModel Editor Demo")
        self.setMinimumSize(700, 500)
        
        # Create LineShapeModel with sample data
        self.lineshape = arts.LineShapeModel()
        self.lineshape.T0 = 296.0  # Reference temperature in Kelvin
        self.lineshape.one_by_one = False
        # single_models is a LineShapeModelMap - starts empty
        
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
        btn = QPushButton("Edit LineShapeModel")
        btn.clicked.connect(self.edit_lineshape)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current values."""
        text = "Current LineShapeModel:\n"
        text += "=" * 60 + "\n\n"
        text += f"Reference Temperature (T0):  {self.lineshape.T0:.1f} K\n"
        text += f"Apply One-by-One:            {self.lineshape.one_by_one}\n"
        text += f"Single Models:               {len(self.lineshape.single_models)} models\n"
        text += "\n" + "=" * 60 + "\n"
        text += "\nLine Shape Model Configuration\n\n"
        text += "Line shape models describe how spectral lines are broadened\n"
        text += "by various physical processes:\n\n"
        text += "Common broadening mechanisms:\n"
        text += "  • Pressure broadening (Lorentz profile)\n"
        text += "  • Doppler broadening (Gaussian profile)\n"
        text += "  • Voigt profile (combination of both)\n"
        text += "  • Dicke narrowing\n"
        text += "  • Speed-dependent effects\n\n"
        text += "T0 is the reference temperature at which line parameters\n"
        text += "are specified (typically 296 K).\n\n"
        text += "Click 'Edit LineShapeModel' to modify configuration.\n"
        text += "Double-click 'single_models' to configure individual\n"
        text += "broadening models for different collision partners."
        self.text.setText(text)
    
    def edit_lineshape(self):
        """Open the LineShapeModel editor."""
        result = edit(self.lineshape, parent=self)
        if result is not None:
            self.lineshape = result
            self.update_display()
            print("\n✓ LineShapeModel updated!")
            print(f"  T0: {self.lineshape.T0:.1f} K")
            print(f"  Models: {len(self.lineshape.single_models)}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = LineShapeModelEditorDemo()
    demo.show()
    sys.exit(app.exec_())
