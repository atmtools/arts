# User-Facing Struct Editors Demo Guide

This guide describes the 5 interactive demos for the newly created user-facing struct editors.

## Demos Overview

All demos follow a consistent pattern:
- Display current values of the struct
- "Edit [Type]" button to open the table-based editor
- Educational descriptions of what each type represents
- Sample data pre-populated for exploration

## 1. ScatteringMetaData Editor (`demo_scattering_metadata_editor.py`)

**Purpose**: Edit metadata for scattering particles used in radiative transfer calculations.

**Fields**:
- `description`: Text description of the particle
- `source`: Data source reference
- `refr_index`: Complex refractive index (string format)
- `mass`: Particle mass in kg
- `diameter_max`: Maximum diameter in meters
- `diameter_volume_equ`: Volume-equivalent diameter in meters
- `diameter_area_equ_aerodynamical`: Aerodynamical area-equivalent diameter in meters

**Sample Data**: Pre-configured ice particle with typical values

**Use Cases**:
- Defining particle properties for cloud/aerosol scattering
- Setting up single scattering databases
- Configuring radar/lidar simulations

---

## 2. ZeemanLineModel Editor (`demo_zeeman_line_model_editor.py`)

**Purpose**: Configure Zeeman effect parameters for spectral lines in magnetic fields.

**Fields**:
- `on`: Enable/disable Zeeman effect (boolean)
- `gu`: Upper state Landé g-factor
- `gl`: Lower state Landé g-factor

**Sample Data**: Pre-configured with typical g-factors (2.002, 1.998)

**Use Cases**:
- Modeling spectral line splitting in Earth's magnetic field
- Solar/stellar atmosphere radiative transfer
- Magnetic field remote sensing

**Background**: The Zeeman effect causes spectral line splitting in magnetic fields. The magnitude of splitting depends on the Landé g-factors of the upper and lower quantum states.

---

## 3. CIARecord Editor (`demo_cia_record_editor.py`)

**Purpose**: Edit Collision-Induced Absorption (CIA) records for molecular pairs.

**Fields**:
- `specs`: Two-element list of SpeciesEnum (read-only after construction)
- `data`: ArrayOfGriddedField2 containing absorption data

**Sample Data**: H2O-N2 pair with empty data array

**Use Cases**:
- Atmospheric absorption modeling (N2-N2, O2-O2)
- Planetary atmosphere studies
- Far-infrared/microwave radiative transfer

**Special Note**: The `specs` field is read-only after construction. Use constructor: `CIARecord(data, species1, species2)`

**Common Pairs**:
- N2-N2: Nitrogen dimer
- H2O-N2: Water-nitrogen complex
- O2-O2: Oxygen dimer
- CO2-CO2: Carbon dioxide dimer

---

## 4. XsecRecord Editor (`demo_xsec_record_editor.py`)

**Purpose**: View cross-section absorption records (read-only structure).

**Fields**:
- `species`: SpeciesEnum
- `version`: Version number
- `fitminpressures`: Vector of minimum pressures
- `fitmaxpressures`: Vector of maximum pressures
- `fitmintemperatures`: Vector of minimum temperatures
- `fitmaxtemperatures`: Vector of maximum temperatures
- `fitcoeffs`: ArrayOfGriddedField4 with polynomial coefficients

**Sample Data**: Empty XsecRecord with default values

**Use Cases**:
- Pre-calculated absorption for complex molecules (CFCs, heavy molecules)
- Fast radiative transfer without line-by-line calculations
- UV-visible band modeling (e.g., O3 Hartley-Huggins bands)

**Special Note**: All XsecRecord properties are read-only after creation. In practice, these objects are loaded from data files rather than constructed manually.

---

## 5. LineShapeModel Editor (`demo_line_shape_model_editor.py`)

**Purpose**: Configure line shape broadening models for spectral lines.

**Fields**:
- `T0`: Reference temperature in Kelvin (typically 296 K)
- `one_by_one`: Apply models one-by-one (boolean)
- `single_models`: LineShapeModelMap with individual broadening models

**Sample Data**: Pre-configured with T0=296 K, one_by_one=False

**Use Cases**:
- Configuring pressure broadening (Lorentz profile)
- Temperature-dependent line widths
- Advanced effects (Dicke narrowing, speed-dependence)
- Collision partner-specific broadening

**Background**: Line shapes describe how spectral lines are broadened by:
- Pressure broadening (collisions → Lorentz profile)
- Doppler broadening (thermal motion → Gaussian profile)
- Combined effects (Voigt profile)
- Advanced effects (Dicke narrowing, speed-dependent broadening)

The `single_models` map allows different broadening parameters for different collision partners (e.g., self-broadening vs. air-broadening).

---

## Running the Demos

All demos can be run directly:

```bash
# Individual demos
python3 demo_scattering_metadata_editor.py
python3 demo_zeeman_line_model_editor.py
python3 demo_cia_record_editor.py
python3 demo_xsec_record_editor.py
python3 demo_line_shape_model_editor.py

# Or run all demos
python3 demo_all_editors.py
```

## Editor Interaction

1. **Launch demo**: Window shows current values with descriptions
2. **Click "Edit" button**: Opens table-based editor dialog
3. **View fields**: Three columns - Field, Type, Value
4. **Edit values**: Double-click any Value cell to edit
5. **Nested editing**: For complex types (arrays, enums), double-clicking dispatches to appropriate sub-editor
6. **Save changes**: Click OK to apply, Cancel to discard

## Technical Notes

### CIARecord Constructor Pattern
Unlike most ARTS types, CIARecord requires constructor arguments:
```python
data = arts.ArrayOfGriddedField2()
cia = arts.CIARecord(data, arts.SpeciesEnum.H2O, arts.SpeciesEnum.N2)
```

The editor still uses the table-based pattern, but be aware that `specs` cannot be modified after construction.

### XsecRecord Read-Only Properties
XsecRecord has all read-only properties. The editor displays the structure but modifications won't persist. In real workflows, XsecRecord objects are loaded from files:
```python
xsec = arts.XsecRecord()
xsec.readxml("xsec_data.xml")
```

### LineShapeModel Advanced Usage
The `single_models` field is a LineShapeModelMap that allows specifying different broadening parameters for each collision partner. This enables accurate modeling of self-broadening vs. foreign-gas broadening.

---

## Coverage Summary

With these 5 editors added:
- **Total editable types**: 256/348 (73.6%)
- **User-facing structs covered**: 5/8 substantive types
- **Remaining**: 3 types (LineShapeModelList and 2 internal array types)

These editors cover the most commonly used scattering and absorption types in ARTS workflows.
