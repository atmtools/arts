# Type Preservation in GUI Editors - Summary

## Problem
When editing ARTS types through the GUI (e.g., PropmatVector, Vector, Matrix, Tensor3), the types were being lost and converted to plain numpy arrays. This broke the ARTS API which expects specific types.

**Root Cause**: The `ndarray.py` editor was wrapping the result in `np.array()`, which stripped away the type that `edit_ndarraylike()` had carefully preserved.

## Solution
Added type preservation logic at multiple dispatch points:

### 1. edit_ndarraylike() - All dimensions (0D, 1D, 2D, 3D+)
**File**: `python/src/pyarts3/gui/common.py`

**Changes**:
- Capture `original_type = type(value)` at the start
- After editing, reconstruct the original type if `original_type != np.ndarray`:
  ```python
  if original_type != np.ndarray:
      try:
          return original_type(result_array)
      except Exception:
          return result_array
  return result_array
  ```

**Applies to**:
- 0D (scalars): Dispatches to `edit.edit()` and reconstructs on return
- 1D arrays: Vector types
- 2D arrays: Matrix, PropmatVector (with constraint: last dim must be 7)
- 3D+ arrays: Tensor3, Tensor4, etc.

### 2. ArrayOf editor - Element editing
**File**: `python/src/pyarts3/gui/edit/ArrayOf.py`

**Changes**:
- Capture `original_type = type(current_obj)` before dispatching
- After `dispatch_edit()` returns, reconstruct if type changed:
  ```python
  if original_type != result_type:
      try:
          new_obj = original_type(new_obj)
      except Exception:
          pass  # Use new_obj as-is
  ```

## How It Works

### Type Preservation Flow
1. **Capture** original type before any conversion
2. **Edit** using numpy arrays for UI convenience
3. **Check** if type changed after editing
4. **Reconstruct** using original type's constructor if needed
5. **Fallback** to edited value if reconstruction fails

### Example: PropmatVector (2D with constraint)
```python
# Original
pv = arts.PropmatVector(np.zeros((3, 7)))  # shape (3, 7)
type(pv)  # PropmatVector

# After editing in GUI
edited_array = np.array([...])  # shape (3, 7), type numpy.ndarray

# Reconstruction
if type(pv) != np.ndarray:
    result = type(pv)(edited_array)  # PropmatVector again!
```

### Example: ArrayOfIndex element (0D scalar)
```python
# Original
element = aoi[1]  # numpy.int64(20)

# Extract for editing
scalar = np.array(element).item()  # int(20)

# After edit.edit() returns
edited = 25  # int

# Reconstruction
if type(element) != type(edited):  # numpy.int64 != int
    result = type(element)(edited)  # numpy.int64(25)
```

## Edge Cases Handled

1. **Type constraints**: PropmatVector requires last dimension = 7
   - Reconstruction fails if shape is wrong
   - Falls back to numpy array (better than crash)

2. **Scalar dispatch**: 0D arrays route through main dispatcher
   - Preserves int vs float distinction
   - Reconstructs numpy.int64, etc.

3. **Failed reconstruction**: Always has fallback
   - Returns edited value even if type conversion fails
   - No data loss

## Verified Test Cases

✅ ArrayOfIndex scalar editing (numpy.int64 preserved)
✅ Vector (1D)
✅ Matrix (2D)
✅ PropmatVector (2D with constraint)
✅ Tensor3 (3D)

## Files Modified

1. `python/src/pyarts3/gui/common.py` - edit_ndarraylike()
   - Added original_type capture
   - Added reconstruction at all return points (1D, 2D, 3D+, 0D)

2. `python/src/pyarts3/gui/edit/ArrayOf.py` - on_cell_double_clicked()
   - Added original_type capture and reconstruction after dispatch

3. `python/src/pyarts3/gui/edit/ndarray.py` - edit()
   - **CRITICAL FIX**: Removed `np.array(result)` wrapper
   - Now returns `edit_ndarraylike(value, parent)` directly
   - Previously was stripping the type that edit_ndarraylike had just preserved!

## Tests

- `test_arrayof_type_preservation.py` - ArrayOfIndex scalar (0D)
- `test_ndarraylike_type_preservation.py` - All dimensions (1D/2D/3D)
