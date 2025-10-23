#!/usr/bin/env python3
"""List all types that do NOT have editors"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from pyarts3 import arts
from pyarts3.utils import builtin_groups
from pyarts3.gui.edit import get_editable_types

print("=" * 80)
print("Finding all NON-EDITABLE types")
print("=" * 80)

# Get all builtin types
all_types = set(builtin_groups())
print(f"\nTotal builtin types: {len(all_types)}")

# Get all editable types
editable_types = get_editable_types()
print(f"Total editable types: {len(editable_types)}")

# Find the difference
non_editable_set = all_types - editable_types
non_editable = sorted([t.__name__ for t in non_editable_set])
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
    elif 'Surface' in t or 'surface' in t.lower():
        categories['Surface'].append(t)
    elif 'Time' in t or 'time' in t.lower():
        categories['Time'].append(t)
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
for category, types in sorted(categories.items()):
    if types:
        print(f"  {category}: {len(types)}")

print("\n" + "=" * 80)
print("Let's check what these types actually are...")
print("=" * 80)

# Sample a few types to understand their structure
sample_types = []
for category, types in categories.items():
    if types:
        sample_types.append((category, types[0]))

for category, typename in sample_types[:10]:  # Check first 10
    try:
        type_obj = getattr(arts, typename)
        print(f"\n{typename} ({category}):")
        
        # Check if it has common attributes
        if hasattr(type_obj, '__dict__'):
            instance = None
            try:
                instance = type_obj()
            except:
                pass
            
            if instance:
                attrs = [a for a in dir(instance) if not a.startswith('_')]
                print(f"  Attributes: {', '.join(attrs[:10])}")
                if len(attrs) > 10:
                    print(f"  ... and {len(attrs) - 10} more")
        
        # Check docstring
        if type_obj.__doc__:
            doc_lines = type_obj.__doc__.strip().split('\n')
            print(f"  Doc: {doc_lines[0][:80]}")
    except Exception as e:
        print(f"\n{typename} ({category}): Error - {e}")
