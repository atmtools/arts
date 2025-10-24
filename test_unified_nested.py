#!/usr/bin/env python3
"""Test unified editor with complex nested types."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication
from pyarts3 import arts
from pyarts3.gui.edit.UnifiedPropertyEditor import edit

app = QApplication(sys.argv)

# Test with AtmPoint (has nested complex types)
print("Testing AtmPoint (has 'mag' which is MagneticAngles)...")
atm = arts.AtmPoint()
atm.temperature = 250.0
atm.pressure = 101325.0

result = edit(atm)

if result is not None:
    print(f"✓ AtmPoint edited!")
    print(f"  temperature: {result.temperature}")
    print(f"  pressure: {result.pressure}")
else:
    print("✗ Cancelled")

sys.exit(0)
