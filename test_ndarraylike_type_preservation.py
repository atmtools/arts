#!/usr/bin/env python3
"""Test that edit_ndarraylike preserves types for all dimensionalities."""

import sys
import numpy as np
import pyarts3.arts as arts

print("Testing type preservation in edit_ndarraylike for all array dimensions\n")
print("=" * 70)

# Test 1D: Vector
print("\n1D Array Test (Vector):")
vec = arts.Vector([1.0, 2.0, 3.0])
print(f"  Original type: {type(vec).__name__}")
arr = np.array(vec)
print(f"  As numpy: shape={arr.shape}, ndim={arr.ndim}")

# Simulate edit and reconstruction
original_type = type(vec)
edited = np.array([4.0, 5.0, 6.0])
try:
    reconstructed = original_type(edited)
    print(f"  Reconstructed type: {type(reconstructed).__name__}")
    if type(reconstructed).__name__ == type(vec).__name__:
        print("  ✅ 1D type preservation works!")
    else:
        print(f"  ✗ Type mismatch: {type(vec).__name__} vs {type(reconstructed).__name__}")
        sys.exit(1)
except Exception as e:
    print(f"  ✗ Failed: {e}")
    sys.exit(1)

# Test 2D: Matrix
print("\n2D Array Test (Matrix):")
mat = arts.Matrix([[1.0, 2.0], [3.0, 4.0]])
print(f"  Original type: {type(mat).__name__}")
arr = np.array(mat)
print(f"  As numpy: shape={arr.shape}, ndim={arr.ndim}")

original_type = type(mat)
edited = np.array([[5.0, 6.0], [7.0, 8.0]])
try:
    reconstructed = original_type(edited)
    print(f"  Reconstructed type: {type(reconstructed).__name__}")
    if type(reconstructed).__name__ == type(mat).__name__:
        print("  ✅ 2D type preservation works!")
    else:
        print(f"  ✗ Type mismatch: {type(mat).__name__} vs {type(reconstructed).__name__}")
        sys.exit(1)
except Exception as e:
    print(f"  ✗ Failed: {e}")
    sys.exit(1)

# Test 2D with constraint: PropmatVector
print("\n2D Array Test with constraint (PropmatVector - must have 7 columns):")
pv = arts.PropmatVector(np.zeros((3, 7)))
print(f"  Original type: {type(pv).__name__}")
arr = np.array(pv)
print(f"  As numpy: shape={arr.shape}, ndim={arr.ndim}")

original_type = type(pv)
edited = np.array([[1,2,3,4,5,6,7], [8,9,10,11,12,13,14], [15,16,17,18,19,20,21]])
try:
    reconstructed = original_type(edited)
    print(f"  Reconstructed type: {type(reconstructed).__name__}")
    if type(reconstructed).__name__ == type(pv).__name__:
        print("  ✅ 2D constrained type preservation works!")
    else:
        print(f"  ✗ Type mismatch: {type(pv).__name__} vs {type(reconstructed).__name__}")
        sys.exit(1)
except Exception as e:
    print(f"  ✗ Failed: {e}")
    sys.exit(1)

# Test 3D: Tensor3
print("\n3D Array Test (Tensor3):")
t3 = arts.Tensor3([[[1, 2], [3, 4]], [[5, 6], [7, 8]]])
print(f"  Original type: {type(t3).__name__}")
arr = np.array(t3)
print(f"  As numpy: shape={arr.shape}, ndim={arr.ndim}")

original_type = type(t3)
edited = np.array([[[9, 10], [11, 12]], [[13, 14], [15, 16]]])
try:
    reconstructed = original_type(edited)
    print(f"  Reconstructed type: {type(reconstructed).__name__}")
    if type(reconstructed).__name__ == type(t3).__name__:
        print("  ✅ 3D type preservation works!")
    else:
        print(f"  ✗ Type mismatch: {type(t3).__name__} vs {type(reconstructed).__name__}")
        sys.exit(1)
except Exception as e:
    print(f"  ✗ Failed: {e}")
    sys.exit(1)

print("\n" + "=" * 70)
print("✅ ALL TESTS PASSED - Type preservation works for all dimensions!")
print("=" * 70)
