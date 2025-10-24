#!/usr/bin/env python3
"""Test the unified property editor with various ARTS types."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication
from pyarts3 import arts
from pyarts3.gui.edit.UnifiedPropertyEditor import edit, inspect_type

# Test inspection
print("Testing inspection...")
print("\n1. AtmPoint:")
atm = arts.AtmPoint()
ro, rw = inspect_type(atm)
print(f"  Read-only: {ro}")
print(f"  Read-write: {rw}")

print("\n2. Sun:")
sun = arts.Sun()
ro, rw = inspect_type(sun)
print(f"  Read-only: {ro}")
print(f"  Read-write: {rw}")

print("\n3. ZeemanLineModel:")
zeeman = arts.ZeemanLineModel()
ro, rw = inspect_type(zeeman)
print(f"  Read-only: {ro}")
print(f"  Read-write: {rw}")

# Test GUI
print("\n\nLaunching GUI test...")
app = QApplication(sys.argv)

# Test with Sun
sun = arts.Sun()
sun.description = "Test sun"
sun.latitude = 45.0
sun.longitude = -122.0

print("\nEditing Sun object...")
result = edit(sun)

if result is not None:
    print(f"\n✓ Edited successfully!")
    print(f"  description: {result.description}")
    print(f"  latitude: {result.latitude}")
    print(f"  longitude: {result.longitude}")
else:
    print("\n✗ Cancelled")

sys.exit(0)
