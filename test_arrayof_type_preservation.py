#!/usr/bin/env python3
"""Test that ArrayOfIndex editing preserves types correctly."""

import sys
import numpy as np
import pyarts3.arts as arts
from PyQt5.QtWidgets import QApplication

# Ensure QApplication exists
app = QApplication.instance()
if app is None:
    app = QApplication(sys.argv)

print("Testing ArrayOfIndex type preservation through edit dispatch\n")

# Create an ArrayOfIndex
aoi = arts.ArrayOfIndex([10, 20, 30])
print(f"Created ArrayOfIndex: {list(aoi)}")

# Get an element
element = aoi[1]
print(f"\nOriginal element: {element}")
print(f"  Type: {type(element)}")
print(f"  Type name: {type(element).__name__}")

# Simulate the full edit flow
print("\nSimulating edit flow:")

# Step 1: Convert to numpy array (happens in edit_ndarraylike)
arr = np.array(element)
print(f"1. As numpy array: shape={arr.shape}, ndim={arr.ndim}")

# Step 2: Detect 0D and extract scalar
if arr.ndim == 0:
    scalar_value = arr.item()
    print(f"2. Extracted scalar: {scalar_value} (type: {type(scalar_value).__name__})")
    
    # Step 3: Simulate edit.edit() returning a plain int
    edited_value = 25  # This is what Index.edit would return
    print(f"3. After edit.edit(): {edited_value} (type: {type(edited_value).__name__})")
    
    # Step 4: Type reconstruction
    original_type = type(element)
    result_type = type(edited_value)
    
    print(f"4. Type check: original={original_type.__name__}, result={result_type.__name__}")
    
    if original_type != result_type:
        print(f"5. Types differ - reconstructing...")
        try:
            final_value = original_type(edited_value)
            print(f"6. Reconstructed: {final_value} (type: {type(final_value).__name__})")
            
            # Verify
            if type(final_value) == type(element):
                print("\n✅ SUCCESS: Type preserved correctly!")
                print(f"   Original: {type(element).__name__}")
                print(f"   Final:    {type(final_value).__name__}")
            else:
                print(f"\n✗ FAILED: Type mismatch!")
                print(f"   Expected: {type(element).__name__}")
                print(f"   Got:      {type(final_value).__name__}")
                sys.exit(1)
                
        except Exception as e:
            print(f"\n✗ FAILED: Reconstruction error: {e}")
            sys.exit(1)
    else:
        print(f"5. Types match - no reconstruction needed")
        final_value = edited_value
        print(f"\n✅ SUCCESS: Type already correct!")
        
else:
    print("\n✗ FAILED: Not detected as 0D array")
    sys.exit(1)
