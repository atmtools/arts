#!/usr/bin/env python3
"""List all types that do NOT have editors - CORRECTED VERSION"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from pyarts3 import arts
from pyarts3.utils import builtin_groups
from pyarts3.gui.edit import get_editable_types

print("=" * 80)
print("Finding all NON-EDITABLE types")
print("=" * 80)

# Get all builtin type names
all_type_names = set(t.__name__ for t in builtin_groups())
print(f"\nTotal builtin types: {len(all_type_names)}")

# Get all editable type names (get_editable_types returns strings already)
editable_type_names = get_editable_types()
print(f"Total editable types: {len(editable_type_names)}")

# Find the difference
non_editable = sorted(all_type_names - editable_type_names)
print(f"Total NON-editable types: {len(non_editable)}")

print("\n" + "=" * 80)
print("NON-EDITABLE TYPES:")
print("=" * 80)

# Try to categorize them
categories = {
    'Agenda': [],
    'Callback': [],
    'Covariance': [],
    'DisortSettings': [],
    'Model': [],
    'Point': [],
    'Quantum': [],
    'Scattering': [],
    'Sensor': [],
    'Spectral': [],
    'Surface': [],
    'Time': [],
    'Workspace': [],
    'Other': []
}

for t in non_editable:
    if 'agenda' in t.lower() or 'Agenda' in t:
        categories['Agenda'].append(t)
    elif 'Callback' in t or 'callback' in t.lower():
        categories['Callback'].append(t)
    elif 'Covariance' in t or 'covariance' in t.lower():
        categories['Covariance'].append(t)
    elif 'Disort' in t or 'disort' in t.lower():
        categories['DisortSettings'].append(t)
    elif 'Model' in t or 'model' in t.lower():
        categories['Model'].append(t)
    elif 'Point' in t or 'point' in t.lower():
        categories['Point'].append(t)
    elif 'Quantum' in t or 'quantum' in t.lower():
        categories['Quantum'].append(t)
    elif 'Scattering' in t or 'scattering' in t.lower() or 'Scat' in t:
        categories['Scattering'].append(t)
    elif 'Sensor' in t or 'sensor' in t.lower():
        categories['Sensor'].append(t)
    elif 'Spectral' in t or 'spectral' in t.lower():
        categories['Spectral'].append(t)
    elif 'Surface' in t or 'surface' in t.lower() or 'Subsurface' in t:
        categories['Surface'].append(t)
    elif 'Time' in t or 'time' in t.lower():
        categories['Time'].append(t)
    elif 'Workspace' in t or 'workspace' in t.lower() or 'Wsv' in t or t in ['Method', 'data', 'parameters', 'line_key']:
        categories['Workspace'].append(t)
    else:
        categories['Other'].append(t)

# Print categorized list
for category, types in sorted(categories.items()):
    if types:
        print(f"\n{category} ({len(types)}):")
        for t in sorted(types):
            print(f"  - {t}")

print("\n" + "=" * 80)
print("SUMMARY BY CATEGORY:")
print("=" * 80)
total = 0
for category, types in sorted(categories.items()):
    if types:
        print(f"  {category}: {len(types)}")
        total += len(types)
print(f"\n  TOTAL: {total}")

print("\n" + "=" * 80)
print("Checking if some can be easily added...")
print("=" * 80)

# Check some specific types that might be easy to add
easy_candidates = []
for typename in non_editable:
    try:
        type_obj = getattr(arts, typename)
        
        # Check if it's a struct-like type with settable attributes
        try:
            instance = type_obj()
            attrs = [a for a in dir(instance) if not a.startswith('_') and not callable(getattr(instance, a))]
            if len(attrs) > 0 and len(attrs) < 20:  # Reasonable number of fields
                # Check if attributes are settable
                settable = []
                for attr in attrs[:5]:  # Check first 5
                    try:
                        val = getattr(instance, attr)
                        setattr(instance, attr, val)
                        settable.append(attr)
                    except:
                        pass
                
                if len(settable) > 0:
                    easy_candidates.append((typename, attrs, settable))
        except:
            pass
    except:
        pass

if easy_candidates:
    print("\nPotential struct-like types that could be easily added:")
    for typename, attrs, settable in easy_candidates[:20]:
        print(f"\n  {typename}:")
        print(f"    All attributes: {', '.join(attrs)}")
        print(f"    Settable: {', '.join(settable)}")
