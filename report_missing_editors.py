#!/usr/bin/env python3
"""Report on missing editors for ARTS types"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from pyarts3 import arts
from pyarts3.utils import builtin_groups
from pyarts3.gui.edit import get_editable_types

# Get all ARTS types
all_types = builtin_groups()
all_type_names = sorted([t.__name__ for t in all_types])

# Get editable types
editable = get_editable_types()

# Find missing types
missing = []
for type_name in all_type_names:
    if type_name not in editable:
        if not type_name.startswith('ArrayOf'):
            missing.append(type_name)

print("=" * 80)
print("ARTS TYPE EDITOR COVERAGE REPORT")
print("=" * 80)
print()
print(f"Total ARTS types:                 {len(all_type_names):3d}")
print(f"Editable types (with editors):    {len(editable):3d} ({100*len(editable)/len(all_type_names):.1f}%)")
print(f"ArrayOf* types (generic editor):  {sum(1 for t in all_type_names if t.startswith('ArrayOf')):3d}")
print(f"Missing editors (non-ArrayOf):    {len(missing):3d}")
print()

# Categorize missing types
print("=" * 80)
print("CATEGORIES OF MISSING TYPES")
print("=" * 80)
print()

categories = {
    'User-Facing Structs (Good candidates for editors)': [
        'Sun', 'PropagationPathPoint', 'SensorObsel', 'SpeciesTag',
        'JacobianTargets', 'DisortSettings', 'AbsorptionBand', 'AbsorptionLine',
        'CIARecord', 'XsecRecord', 'ScatteringMetaData', 'ZeemanLineModel',
        'LineShapeModel'
    ],
    'Internal/Advanced Types (Low priority)': [
        'Agenda', 'Method', 'Block', 'CallbackOperator', 'Workspace*',
        'disort_settings_agendaOperator', 'propagation_matrix_agendaOperator',
        'spectral_radiance_*_agendaOperator', 'surface_reflectance_agendaOperator',
        'CxxWorkspace', 'Wsv'
    ],
    'Operator Types (Usually auto-generated)': [
        'NumericBinaryOperator', 'NumericTernaryOperator', 'NumericUnaryOperator',
        'SpectralRadianceOperator', 'SpectralRadianceTransformOperator',
        'ExtSSACallback', 'DisortBDRFOperator'
    ],
    'Grid Types (Already have array-like editors)': [
        'DoubleGaussGrid', 'FejerGrid', 'GaussLegendreGrid', 'LobattoGrid',
        'IrregularZenithAngleGrid'
    ],
    'Complex Scattering Types': [
        'SingleScatteringData*', 'BulkScatteringProperties*', 
        'ScatteringGeneralSpectral*', 'HenyeyGreensteinScatterer',
        'ParticleHabit', 'ParticleProperties', 'BinnedPSD'
    ],
    'Specialized Data Structures': [
        'AbsorptionLookupTable', 'CovarianceMatrix', 'PartitionFunctionsData',
        'MTCKD400WaterData', 'SpeciesIsotopologueRatios', 'Sparse',
        'MatrixOfDisortBDRF', 'DisortFlux', 'DisortRadiance'
    ],
    'Tag/Enum Types (Usually edited via Options editor)': [
        'SubsurfacePropertyTag', 'SurfacePropertyTag', 'ErrorKey',
        'JacobianTargetType', 'PolarizationState', 'PType', 'ModelName'
    ]
}

for category, types in categories.items():
    matching = [t for t in missing if any(pattern.replace('*', '') in t for pattern in types)]
    if matching:
        print(f"{category}:")
        print(f"  Count: {len(matching)}")
        print(f"  Examples: {', '.join(matching[:5])}")
        if len(matching) > 5:
            print(f"           ... and {len(matching) - 5} more")
        print()

print("=" * 80)
print("RECOMMENDED PRIORITIES FOR NEW EDITORS")
print("=" * 80)
print()

high_priority = [
    ('Sun', 'Solar properties - distance, position, spectrum'),
    ('PropagationPathPoint', 'Ray path point with position and line-of-sight'),
    ('SpeciesTag', 'Species identification with type and CIA info'),
]

medium_priority = [
    ('SensorObsel', 'Sensor observation element'),
    ('JacobianTargets', 'Jacobian calculation target configuration'),
    ('DisortSettings', 'DISORT radiative transfer settings'),
    ('AbsorptionBand', 'Absorption band with lines and lineshape'),
    ('AbsorptionLine', 'Individual absorption line parameters'),
]

print("HIGH PRIORITY (Simple user-facing types):")
for i, (type_name, desc) in enumerate(high_priority, 1):
    status = "✓ Has editor" if type_name in editable else "✗ Missing"
    print(f"  {i}. {type_name:25s} - {desc:45s} [{status}]")

print()
print("MEDIUM PRIORITY (More complex but user-facing):")
for i, (type_name, desc) in enumerate(medium_priority, 1):
    status = "✓ Has editor" if type_name in editable else "✗ Missing"
    print(f"  {i}. {type_name:25s} - {desc:45s} [{status}]")

print()
print("=" * 80)
print(f"SUMMARY: {len(editable)} of {len(all_type_names)} types have editors ({100*len(editable)/len(all_type_names):.1f}%)")
print("=" * 80)
