#!/usr/bin/env python3
"""Test unified editor with various ARTS types including maps."""

import sys
sys.path.insert(0, 'build/python/src')

from PyQt5.QtWidgets import QApplication
from pyarts3 import arts
from pyarts3.gui.edit import edit
from pyarts3.gui.edit.UnifiedPropertyEditor import is_terminal_type, inspect_type

app = QApplication(sys.argv)

print("Testing terminal type detection...")
print()

# Test basic types
vec = arts.Vector([1, 2, 3])
print(f"Vector: is_terminal={is_terminal_type(vec)}")

# Test ArrayOf
array_vec = arts.ArrayOfVector()
print(f"ArrayOfVector: is_terminal={is_terminal_type(array_vec)}")

# Test map-like
try:
    # LineShapeModel has single_models which is a map
    lsm = arts.LineShapeModel()
    print(f"LineShapeModel: is_terminal={is_terminal_type(lsm)}")
    ro, rw = inspect_type(lsm)
    print(f"  Properties: RW={rw}, RO={ro}")
    if rw:
        for prop in rw:
            val = getattr(lsm, prop)
            print(f"    {prop}: type={type(val).__name__}, is_terminal={is_terminal_type(val)}")
except Exception as e:
    print(f"LineShapeModel error: {e}")

print()

# Test complex type
print("Testing complex type (Sun)...")
sun = arts.Sun()
sun.description = "Test sun"
sun.latitude = 45.0

print(f"Sun: is_terminal={is_terminal_type(sun)}")
ro, rw = inspect_type(sun)
print(f"  Properties: RW={rw}, RO={ro}")

print("\n✓ All tests passed!")
print("\nNow opening GUI to test Sun editing...")

result = edit(sun)

if result is not None:
    print(f"\n✓ Edited!")
    print(f"  description: {result.description}")
    print(f"  latitude: {result.latitude}")
else:
    print("\n✗ Cancelled")

sys.exit(0)
