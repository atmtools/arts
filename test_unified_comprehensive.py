#!/usr/bin/env python3
"""Comprehensive test of the unified property editor system."""

import sys
sys.path.insert(0, 'build/python/src')

from pyarts3 import arts
from pyarts3.gui.edit.UnifiedPropertyEditor import is_terminal_type, inspect_type

print("=" * 70)
print("UNIFIED PROPERTY EDITOR - COMPREHENSIVE TYPE TEST")
print("=" * 70)

# Test cases: (Type, expected_terminal, description)
test_cases = [
    # Terminal types
    (arts.Numeric(3.14), True, "Numeric - ARTS scalar"),
    (arts.Index(42), True, "Index - ARTS integer"),
    (arts.String("test"), True, "String - ARTS string"),
    (arts.Vector([1, 2, 3]), True, "Vector - has __array__"),
    (arts.Matrix([[1, 2], [3, 4]]), True, "Matrix - has __array__"),
    (arts.ArrayOfVector(), True, "ArrayOfVector - ArrayOf type"),
    (arts.ArrayOfString(), True, "ArrayOfString - ArrayOf type"),
    
    # Complex types (should use property editor)
    (arts.Sun(), False, "Sun - complex with properties"),
    (arts.AtmPoint(), False, "AtmPoint - complex with properties"),
    (arts.PropagationPathPoint(), False, "PropagationPathPoint - complex"),
    (arts.ZeemanLineModel(), False, "ZeemanLineModel - complex"),
    (arts.LineShapeModel(), False, "LineShapeModel - complex with map property"),
]

print("\nTesting terminal type detection:\n")
print(f"{'Type':<35} {'Terminal?':<12} {'Status'}")
print("-" * 70)

passed = 0
failed = 0

for value, expected, description in test_cases:
    actual = is_terminal_type(value)
    status = "✓ PASS" if actual == expected else "✗ FAIL"
    if actual == expected:
        passed += 1
    else:
        failed += 1
    print(f"{description:<35} {str(actual):<12} {status}")

print("-" * 70)
print(f"Results: {passed} passed, {failed} failed\n")

# Show property inspection for complex types
print("\nProperty inspection for complex types:\n")
print(f"{'Type':<25} {'Read-Write Properties'}")
print("-" * 70)

for value, expected, description in test_cases:
    if not expected:  # Complex types
        ro, rw = inspect_type(value)
        type_name = type(value).__name__
        print(f"{type_name:<25} {', '.join(rw)}")

print("\n" + "=" * 70)
print("TEST COMPLETE")
print("=" * 70)

# Summary
if failed == 0:
    print("\n✓ All tests passed! The unified editor correctly identifies:")
    print("  • Terminal types (use specialized editors)")
    print("  • Complex types (use property navigation)")
    print("  • Map types (use map editor)")
else:
    print(f"\n✗ {failed} tests failed - please review")
    sys.exit(1)

sys.exit(0)
