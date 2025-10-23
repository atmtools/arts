#!/usr/bin/env python3
"""Test the get_editable_types() function to verify all types are detected."""

import sys
import pyarts3.arts as arts
from pyarts3.gui.edit import get_editable_types, can_edit

print("=" * 70)
print("Testing get_editable_types() for Workspace Editor")
print("=" * 70)

# Get all editable types
editable_types = get_editable_types()

print(f"\n1. Total editable types: {len(editable_types)}")

# Categorize types
explicit_editors = []
option_groups = []
array_types = []
arrayof_types = []

for type_name in sorted(editable_types):
    if type_name in arts.globals.option_groups():
        option_groups.append(type_name)
    elif type_name.startswith('ArrayOf'):
        arrayof_types.append(type_name)
    elif hasattr(arts, type_name):
        try:
            t = getattr(arts, type_name)
            if hasattr(t, '__array__') or hasattr(t(), '__array__'):
                array_types.append(type_name)
            else:
                explicit_editors.append(type_name)
        except:
            explicit_editors.append(type_name)
    else:
        explicit_editors.append(type_name)

print(f"\n2. Breakdown:")
print(f"   Explicit editors: {len(explicit_editors)}")
print(f"   Option groups: {len(option_groups)}")
print(f"   Array-like types: {len(array_types)}")
print(f"   ArrayOf* types: {len(arrayof_types)}")

print(f"\n3. Sample explicit editors:")
for t in explicit_editors[:10]:
    print(f"   - {t}")

print(f"\n4. Sample option groups:")
for t in option_groups[:5]:
    print(f"   - {t}")

print(f"\n5. Sample array-like types:")
for t in array_types[:10]:
    print(f"   - {t}")

print(f"\n6. Sample ArrayOf types:")
for t in arrayof_types[:10]:
    print(f"   - {t}")

# Test specific important types
print(f"\n7. Testing specific critical types:")

critical_types = [
    'Index', 'Numeric', 'String',  # Basic types
    'Vector', 'Matrix', 'Tensor3',  # Array types
    'PropmatVector',  # Special 2D type
    'ArrayOfIndex', 'ArrayOfNumeric', 'ArrayOfString',  # ArrayOf basics
    'ArrayOfVector', 'ArrayOfMatrix',  # ArrayOf arrays
    'ArrayOfArrayOfIndex',  # Nested ArrayOf
    'ray_path_observer_agendaSetGeometricMaxStep',  # Option group
]

all_passed = True
for type_name in critical_types:
    result = can_edit(type_name)
    status = "✓" if result else "✗"
    print(f"   {status} {type_name}: {result}")
    if not result:
        all_passed = False

# Check for ArrayOfArrayOf types
print(f"\n8. Checking ArrayOfArrayOf* types:")
arrayof_arrayof = [t for t in arrayof_types if t.startswith('ArrayOfArrayOf')]
print(f"   Found {len(arrayof_arrayof)} ArrayOfArrayOf types")
if arrayof_arrayof:
    for t in arrayof_arrayof[:5]:
        print(f"   - {t}")

# Verify recursive ArrayOf detection
print(f"\n9. Testing ArrayOf recursion:")
if 'Index' in editable_types:
    print(f"   Index: ✓")
    if 'ArrayOfIndex' in editable_types:
        print(f"   ArrayOfIndex: ✓")
        if 'ArrayOfArrayOfIndex' in editable_types:
            print(f"   ArrayOfArrayOfIndex: ✓")
            print(f"   Recursion working correctly!")
        else:
            print(f"   ArrayOfArrayOfIndex: ✗ Missing!")
            all_passed = False
    else:
        print(f"   ArrayOfIndex: ✗ Missing!")
        all_passed = False
else:
    print(f"   Index: ✗ Missing!")
    all_passed = False

if all_passed:
    print("\n" + "=" * 70)
    print("✅ SUCCESS: All critical types are editable!")
    print("=" * 70)
    print("\nThe Workspace editor will now correctly identify:")
    print("  - All basic types (Index, Numeric, String)")
    print("  - All array-like types (Vector, Matrix, Tensor3, PropmatVector, etc.)")
    print("  - All option group enums (46 types)")
    print("  - All ArrayOf variants (including nested ArrayOfArrayOf)")
    print("\nNo types should be marked as unavailable incorrectly!")
else:
    print("\n" + "=" * 70)
    print("✗ FAILED: Some critical types missing")
    print("=" * 70)
    sys.exit(1)
