# Demo Scripts for New Editors

## Quick Start

Run the master launcher to select any demo interactively:
```bash
python3 demo_all_editors.py
```

Or run individual demos directly:

## Individual Demo Scripts

### 1. Sun Editor (`demo_sun_editor.py`)
```bash
python3 demo_sun_editor.py
```
**Shows**: Solar properties configuration
- Distance, position (lat/lon), radius
- Solar spectrum (Matrix)
- Interactive editing of all 6 fields

### 2. PropagationPathPoint Editor (`demo_propagation_path_point_editor.py`)
```bash
python3 demo_propagation_path_point_editor.py
```
**Shows**: Ray path point with position and line-of-sight
- 3D position vector (Vector3)
- Line-of-sight angles (Vector2)
- Refractive indices (nreal, ngroup)
- Position/LOS type enums

### 3. SpeciesTag Editor (`demo_species_tag_editor.py`)
```bash
python3 demo_species_tag_editor.py
```
**Shows**: Species identification
- Species enum (H2O, O2, N2, etc.)
- Tag type (Plain, Predefined, CIA, etc.)
- CIA 2nd species
- Full name string

### 4. SensorObsel Editor (`demo_sensor_obsel_editor.py`)
```bash
python3 demo_sensor_obsel_editor.py
```
**Shows**: Sensor observation element
- Frequency grid (AscendingGrid)
- Position/LOS vector (SensorPosLosVector)
- Weight matrix (SparseStokvecMatrix)

### 5. JacobianTargets Editor (`demo_jacobian_targets_editor.py`)
```bash
python3 demo_jacobian_targets_editor.py
```
**Shows**: Jacobian calculation targets
- 6 target arrays (atm, surf, subsurf, line, sensor, error)
- Each is an ArrayOfJacobianTarget
- Used for retrieval/sensitivity analysis

### 6. DisortSettings Editor (`demo_disort_settings_editor.py`)
```bash
python3 demo_disort_settings_editor.py
```
**Shows**: DISORT radiative transfer settings (most complex - 16 fields!)
- Dimensions (quadrature, Fourier, Legendre)
- Altitude and frequency grids
- Optical properties (thickness, albedo, scattering)
- Legendre coefficients (Tensor3)
- Solar parameters (source, angles)
- Boundary conditions (Tensor3)
- BRDF matrix

### 7. AbsorptionBand Editor (`demo_absorption_band_editor.py`)
```bash
python3 demo_absorption_band_editor.py
```
**Shows**: Absorption band with lines
- Cutoff type and value
- Array of absorption lines
- Line shape configuration

### 8. AbsorptionLine Editor (`demo_absorption_line_editor.py`)
```bash
python3 demo_absorption_line_editor.py
```
**Shows**: Individual spectroscopic line
- Center frequency (f0), strength (a), energy (e0)
- Degeneracies (gu, gl)
- Line shape model (ls)
- Zeeman model (z)
- Quantum state (qn)
- Example: H2O line at 183.31 GHz

## Features Demonstrated

All demos show:
- ✅ **Table-based UI** with Field/Type/Value columns
- ✅ **Double-click editing** - click Value to open nested editors
- ✅ **Recursive editing** - edit complex nested types (vectors, matrices, arrays)
- ✅ **OK/Cancel workflow** - commit or rollback changes
- ✅ **Live updates** - see changes reflected immediately

## What to Try

1. **Launch a demo**: Run any script above
2. **Click "Edit" button**: Opens the table-based editor dialog
3. **Double-click a Value cell**: Opens appropriate sub-editor
   - Numeric values → simple input
   - Enums → dropdown selection
   - Vectors/Matrices → array editors
   - Nested structs → recursive struct editors
4. **Click OK**: Commits changes and updates display
5. **Click Cancel**: Discards changes

## Tips

- **Complex types**: DisortSettings shows the most fields (16)
- **Nested editing**: Try editing Vector3 in PropagationPathPoint
- **Arrays**: JacobianTargets shows multiple array fields
- **Enums**: SpeciesTag demonstrates enum selection
- **Scientific data**: AbsorptionLine shows real spectroscopic parameters

## Next Steps

After trying the demos:
1. Use these patterns in your own ARTS workflows
2. Check coverage with: `python3 report_missing_editors.py`
3. See `PRIORITY_EDITORS_SUMMARY.md` for implementation details
