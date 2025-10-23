#!/usr/bin/env python3
"""Interactive demo for the DisortSettings editor."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTextEdit
from pyarts3 import arts
from pyarts3.gui.edit import edit


class DisortSettingsEditorDemo(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DisortSettings Editor Demo")
        self.setMinimumSize(750, 600)
        
        # Create DisortSettings with sample configuration
        self.settings = arts.DisortSettings()
        self.settings.quadrature_dimension = 16
        self.settings.fourier_mode_dimension = 4
        self.settings.legendre_polynomial_dimension = 8
        
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
        btn = QPushButton("Edit DisortSettings")
        btn.clicked.connect(self.edit_settings)
        layout.addWidget(btn)
    
    def update_display(self):
        """Update the text display with current values."""
        text = "Current DisortSettings:\n"
        text += "=" * 50 + "\n\n"
        text += "DISORT Radiative Transfer Configuration:\n\n"
        text += f"Dimensions:\n"
        text += f"  Quadrature:   {self.settings.quadrature_dimension}\n"
        text += f"  Fourier mode: {self.settings.fourier_mode_dimension}\n"
        text += f"  Legendre:     {self.settings.legendre_polynomial_dimension}\n\n"
        text += f"Grids:\n"
        text += f"  Altitude:     {len(self.settings.altitude_grid)} levels\n"
        text += f"  Frequency:    {len(self.settings.frequency_grid)} points\n\n"
        text += f"Optical Properties:\n"
        text += f"  Optical thickness:    {self.settings.optical_thicknesses.shape}\n"
        text += f"  Single scat. albedo:  {self.settings.single_scattering_albedo.shape}\n"
        text += f"  Fractional scatting:  {self.settings.fractional_scattering.shape}\n\n"
        text += f"Solar:\n"
        text += f"  Source:       {len(self.settings.solar_source)} elements\n"
        text += f"  Zenith angle: {len(self.settings.solar_zenith_angle)} elements\n"
        text += f"  Azimuth:      {len(self.settings.solar_azimuth_angle)} elements\n"
        text += "\n" + "=" * 50 + "\n"
        text += "\nClick 'Edit DisortSettings' to open the editor.\n"
        text += "This is the most complex editor with 16 configurable fields!"
        self.text.setText(text)
    
    def edit_settings(self):
        """Open the DisortSettings editor."""
        result = edit(self.settings, parent=self)
        if result is not None:
            self.settings = result
            self.update_display()
            print("\n✓ DisortSettings updated!")
            print(f"  Quadrature dimension: {self.settings.quadrature_dimension}")
            print(f"  Fourier modes: {self.settings.fourier_mode_dimension}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = DisortSettingsEditorDemo()
    demo.show()
    sys.exit(app.exec_())
