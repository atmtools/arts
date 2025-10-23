# Empty Array Handling Fix

## Problem

The DisortSettings demo (and other editors) crashed when double-clicking on empty multi-dimensional arrays (Tensor3, Matrix, etc.) because the reshape operation tried to reshape a 0-sized array:

```
ValueError: cannot reshape array of size 0 into shape (0)
```

This occurred because:
1. Empty Tensor3 has shape `(0, 0, 0)` with size 0
2. Empty Matrix has shape `(0, 0)` with size 0
3. The `edit_ndarraylike` function tried to reshape these without checking if they were empty

## Solution

Added empty array checks in `python/src/pyarts3/gui/common.py`:

### For 2D Arrays (Matrix)
```python
# Handle empty arrays
if data.size == 0:
    layout.addWidget(QLabel(f"Empty array with shape: {data.shape}"))
    layout.addWidget(QLabel("Cannot edit empty arrays."))
    layout.addWidget(QLabel("Please reshape or populate the array first."))
    buttons = QDialogButtonBox(QDialogButtonBox.Ok)
    buttons.accepted.connect(dialog.accept)
    layout.addWidget(buttons)
    dialog.exec_()
    return None
```

### For 3D+ Arrays (Tensor3, etc.)
Same check added before attempting reshape operation.

## Result

- Empty arrays now show a friendly message instead of crashing
- User is informed that the array is empty and cannot be edited until populated
- Non-empty arrays continue to work normally
- All demos (especially DisortSettings) now work without crashes

## Testing

Run the test script:
```bash
python3 test_empty_array_handling.py
```

Or test interactively:
```bash
python3 demo_disort_settings_editor.py
# Click "Edit DisortSettings"
# Double-click on any empty Tensor3 field (like legendre_coefficients)
# Should see friendly message instead of crash
```

## Files Modified

- `python/src/pyarts3/gui/common.py` - Added empty array checks
- `build/python/src/pyarts3/gui/common.py` - Synced

## Impact

This fix prevents crashes when users interact with any struct type that has empty multi-dimensional array fields, including:
- DisortSettings (multiple Tensor3 and Matrix fields)
- SensorObsel (SparseStokvecMatrix)
- Any other type with Matrix/Tensor fields
