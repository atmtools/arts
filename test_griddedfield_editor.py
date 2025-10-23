#!/usr/bin/env python3
"""Test GriddedField editor"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from pyarts3 import arts
import numpy as np

# Test with GriddedField2
print("Testing GriddedField2 editor...")
gf = arts.GriddedField2()
gf.grids = [arts.Vector([1.0, 2.0, 3.0]), arts.Vector([10.0, 20.0])]
gf.gridnames = ['frequency', 'temperature']
gf.dataname = 'absorption'
gf.data = arts.Matrix([[1.1, 1.2], [2.1, 2.2], [3.1, 3.2]])

print(f"Original GriddedField2:")
print(f"  Type: {type(gf).__name__}")
print(f"  Dim: {arts.GriddedField2.dim}")
print(f"  Grids: {[list(g) for g in gf.grids]}")
print(f"  Grid names: {gf.gridnames}")
print(f"  Data name: {gf.dataname}")
print(f"  Data shape: {gf.data.shape}")
print(f"  Data:\n{gf.data}")

# Check detection
print("\nChecking gridded field detection:")
print(f"  hasattr(__array__): {hasattr(gf, '__array__')}")
print(f"  hasattr(grids): {hasattr(gf, 'grids')}")
print(f"  hasattr(gridnames): {hasattr(gf, 'gridnames')}")
print(f"  hasattr(dataname): {hasattr(gf, 'dataname')}")
print(f"  hasattr(data): {hasattr(gf, 'data')}")

# This will open the GUI with tabs for Metadata, Grids, and Data
from pyarts3.gui.edit import edit
result = edit(gf)

if result is not None:
    print(f"\nResult type: {type(result).__name__}")
    print(f"  Grids: {[list(g) for g in result.grids]}")
    print(f"  Grid names: {result.gridnames}")
    print(f"  Data name: {result.dataname}")
    print(f"  Data shape: {result.data.shape}")
    print(f"  Type preserved: {type(result) == type(gf)}")
else:
    print("\nCancelled")

# Test with GriddedField3
print("\n" + "="*60)
print("Testing GriddedField3 editor...")
gf3 = arts.GriddedField3()
gf3.grids = [
    arts.Vector([1.0, 2.0]),
    arts.Vector([10.0, 20.0, 30.0]),
    arts.Vector([100.0, 200.0])
]
gf3.gridnames = ['freq', 'temp', 'pressure']
gf3.dataname = 'cross_section'
gf3.data = arts.Tensor3([
    [[1, 2], [3, 4], [5, 6]],
    [[7, 8], [9, 10], [11, 12]]
])

print(f"Original GriddedField3:")
print(f"  Type: {type(gf3).__name__}")
print(f"  Dim: {arts.GriddedField3.dim}")
print(f"  Grids: {[list(g) for g in gf3.grids]}")
print(f"  Grid names: {gf3.gridnames}")
print(f"  Data name: {gf3.dataname}")
print(f"  Data shape: {gf3.data.shape}")

result3 = edit(gf3)

if result3 is not None:
    print(f"\nResult type: {type(result3).__name__}")
    print(f"  Grid names: {result3.gridnames}")
    print(f"  Data name: {result3.dataname}")
    print(f"  Type preserved: {type(result3) == type(gf3)}")
else:
    print("\nCancelled")
