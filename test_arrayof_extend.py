#!/usr/bin/env python3
"""Test ArrayOf editor with add/remove functionality"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from pyarts3 import arts

# Test with ArrayOfIndex
print("Testing ArrayOfIndex with add/remove...")
arr = arts.ArrayOfIndex([1, 2, 3])
print(f"Original: {list(arr)}")

# This will open the GUI - you should be able to:
# 1. Double-click to edit elements
# 2. Click "Add Element" to add new elements
# 3. Click "Remove Selected" to delete elements
from pyarts3.gui.edit.ArrayOf import edit
result = edit(arr)

if result is not None:
    print(f"Result type: {type(result).__name__}")
    print(f"Result: {list(result)}")
    print(f"Type preserved: {type(result) == type(arr)}")
else:
    print("Cancelled")

print("\nTesting ArrayOfString with add/remove...")
arr2 = arts.ArrayOfString(["hello", "world"])
print(f"Original: {list(arr2)}")

result2 = edit(arr2)

if result2 is not None:
    print(f"Result type: {type(result2).__name__}")
    print(f"Result: {list(result2)}")
    print(f"Type preserved: {type(result2) == type(arr2)}")
else:
    print("Cancelled")
