# CIARecord and XsecRecord Editor Fixes

## Summary

Fixed the CIARecord, XsecRecord, and ZeemanLineModel editors to properly respect the C++ interface definitions in `py_cia.cpp` and `py_xsec_fit.cpp`.

## Root Causes

### 1. CIARecord Interface (py_cia.cpp)

```cpp
.def_prop_ro("specs", ...)  // READ-ONLY property
.def_prop_ro("data", ...)   // READ-ONLY property
```

**Problem**: Editor was trying to use `setattr(result, "specs", ...)` and `setattr(result, "data", ...)` on read-only properties.

**Solution**: Use constructor `CIARecord(data, species1, species2)` to create object with values.

### 2. XsecRecord Interface (py_xsec_fit.cpp)

```cpp
.def_ro_static("version", ...)  // READ-ONLY STATIC property
.def_rw("species", ...)         // READ-WRITE property
.def_rw("fitcoeffs", ...)       // READ-WRITE property
... (all other fields are def_rw)
```

**Problem**: Editor was trying to set `version` which is a read-only class static property.

**Solution**: Skip the `version` field when setting attributes (it's always 2).

### 3. ZeemanLineModel Boolean Type Issue

**Problem**: The `on` field expects a `bool` but Generic editor returns `int` for boolean values.

**Solution**: Explicitly convert to `bool` before setting: `val = bool(val)`.

## Code Changes

### CIARecord.py
```python
# OLD - Tried to setattr on read-only properties
if dialog.exec_() == QDialog.Accepted:
    result = arts.CIARecord()
    for field_name, _, _ in fields:
        setattr(result, field_name, current_values[field_name])
    return result

# NEW - Use constructor
if dialog.exec_() == QDialog.Accepted:
    try:
        specs = current_values["specs"]
        data = current_values["data"]
        result = arts.CIARecord(data, specs[0], specs[1])
        return result
    except Exception as e:
        QMessageBox.warning(dialog, "Error", f"Failed to create CIARecord: {e}")
        return None
```

### XsecRecord.py
```python
# OLD - Tried to set read-only static version
if dialog.exec_() == QDialog.Accepted:
    result = arts.XsecRecord()
    for field_name, _, _ in fields:
        setattr(result, field_name, current_values[field_name])
    return result

# NEW - Skip version field
if dialog.exec_() == QDialog.Accepted:
    result = arts.XsecRecord()
    for field_name, _, _ in fields:
        if field_name != "version":  # Skip read-only version field
            setattr(result, field_name, current_values[field_name])
    return result
```

### ZeemanLineModel.py
```python
# OLD - int passed to bool parameter
if dialog.exec_() == QDialog.Accepted:
    result = arts.ZeemanLineModel()
    for field_name, _, _ in fields:
        setattr(result, field_name, current_values[field_name])
    return result

# NEW - Explicit bool conversion
if dialog.exec_() == QDialog.Accepted:
    result = arts.ZeemanLineModel()
    for field_name, _, _ in fields:
        val = current_values[field_name]
        if field_name == "on":
            val = bool(val)  # Generic editor returns int for bool
        setattr(result, field_name, val)
    return result
```

## C++ Interface Reference

### CIARecord Constructor
```cpp
.def(py::init<ArrayOfGriddedField2, SpeciesEnum, SpeciesEnum>())
```

Creates a CIARecord with:
1. `data`: ArrayOfGriddedField2
2. `species1`: First SpeciesEnum
3. `species2`: Second SpeciesEnum

The `specs` property returns a tuple of the two species.

### XsecRecord Properties

- **version**: Read-only class static (always 2)
- **species**: Read-write SpeciesEnum
- **fitminpressures**: Read-write Vector
- **fitmaxpressures**: Read-write Vector
- **fitmintemperatures**: Read-write Vector
- **fitmaxtemperatures**: Read-write Vector  
- **fitcoeffs**: Read-write ArrayOfGriddedField1Named

## Testing

All 5 user-facing struct demos now work correctly:

```bash
python3 demo_scattering_metadata_editor.py  # ✓ Works
python3 demo_zeeman_line_model_editor.py    # ✓ Works (bool fixed)
python3 demo_cia_record_editor.py           # ✓ Works (constructor)
python3 demo_xsec_record_editor.py          # ✓ Works (skip version)
python3 demo_line_shape_model_editor.py     # ✓ Works
```

## Lessons Learned

1. **Always check C++ binding files** (py_*.cpp) to understand:
   - Which properties are read-only (`def_prop_ro`, `def_ro_static`)
   - Which properties are read-write (`def_rw`, `def_prop_rw`)
   - What constructor signatures are available (`def(py::init<...>())`)

2. **Property patterns**:
   - `def_rw`: Direct read-write access to member variable
   - `def_prop_rw`: Custom getter/setter, usually writable
   - `def_prop_ro`: Custom getter only, READ-ONLY
   - `def_ro_static`: Read-only class static property

3. **Generic editor limitations**:
   - Returns `int` for boolean values (needs explicit `bool()` conversion)
   - Should be enhanced to preserve bool type in future

4. **Constructor-based types**:
   - Some types (like CIARecord) require constructor initialization
   - Cannot modify properties after construction
   - Editor must collect all values and reconstruct object

## Files Modified

- `python/src/pyarts3/gui/edit/CIARecord.py` - Use constructor pattern
- `python/src/pyarts3/gui/edit/XsecRecord.py` - Skip read-only version field
- `python/src/pyarts3/gui/edit/ZeemanLineModel.py` - Convert bool types
- `build/python/src/pyarts3/gui/edit/*.py` - Synced all changes

All editors now properly respect the C++ interface contracts.
