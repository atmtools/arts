#!/usr/bin/env python3
"""Test map-like editor"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from pyarts3 import arts

# Test with QuantumIdentifierNumericMap
print("Testing QuantumIdentifierNumericMap editor...")
m = arts.QuantumIdentifierNumericMap()

# Add some entries
qi1 = arts.QuantumIdentifier('O2-66')
qi2 = arts.QuantumIdentifier('H2O')
qi3 = arts.QuantumIdentifier('CO2')

m[qi1] = 1.5
m[qi2] = 2.5
m[qi3] = 3.7

print(f"Original map:")
print(f"  Type: {type(m).__name__}")
print(f"  Length: {len(m)}")
for k, v in m.items():
    print(f"    {k} -> {v}")

# Check detection
print("\nChecking map-like detection:")
print(f"  hasattr(keys): {hasattr(m, 'keys')}")
print(f"  hasattr(items): {hasattr(m, 'items')}")
print(f"  hasattr(values): {hasattr(m, 'values')}")
print(f"  hasattr(__getitem__): {hasattr(m, '__getitem__')}")
print(f"  hasattr(__setitem__): {hasattr(m, '__setitem__')}")

# This will open the GUI with table showing key-value pairs
# Double-click to edit keys or values
# Use + Add / − Remove buttons to add/remove entries
from pyarts3.gui.edit import edit
result = edit(m)

if result is not None:
    print(f"\nResult type: {type(result).__name__}")
    print(f"  Length: {len(result)}")
    for k, v in result.items():
        print(f"    {k} -> {v}")
    print(f"  Type preserved: {type(result) == type(m)}")
else:
    print("\nCancelled")

# Test with a simpler map type
print("\n" + "="*60)
print("Testing AbsorptionBands (another map-like type)...")
bands = arts.AbsorptionBands()

# Add some bands
from pyarts3.arts import AbsorptionBand, QuantumIdentifier
band1 = arts.AbsorptionBand()
band2 = arts.AbsorptionBand()

qi_a = arts.QuantumIdentifier('O2-66')
qi_b = arts.QuantumIdentifier('H2O')

bands[qi_a] = band1
bands[qi_b] = band2

print(f"Original AbsorptionBands:")
print(f"  Type: {type(bands).__name__}")
print(f"  Length: {len(bands)}")
print(f"  Keys: {list(bands.keys())}")

result2 = edit(bands)

if result2 is not None:
    print(f"\nResult type: {type(result2).__name__}")
    print(f"  Length: {len(result2)}")
    print(f"  Type preserved: {type(result2) == type(bands)}")
else:
    print("\nCancelled")
