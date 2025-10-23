#!/usr/bin/env python3
"""Test QuantumIdentifier editor"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from pyarts3 import arts
from pyarts3.gui.edit import edit, get_editable_types

# Check if QuantumIdentifier is now editable
editable = get_editable_types()
print(f"QuantumIdentifier editable: {'QuantumIdentifier' in editable}")

# Test with QuantumIdentifier
print("\nTesting QuantumIdentifier editor...")
qi = arts.QuantumIdentifier('O2-66 J 1 2')

print(f"Original QuantumIdentifier:")
print(f"  String: {qi}")
print(f"  isot: {qi.isot}")
print(f"  state: {qi.state}")

# This will open a simple dialog to edit the string representation
result = edit(qi)

if result is not None:
    print(f"\nResult:")
    print(f"  String: {result}")
    print(f"  isot: {result.isot}")
    print(f"  state: {result.state}")
    print(f"  Type preserved: {type(result) == type(qi)}")
else:
    print("\nCancelled")
