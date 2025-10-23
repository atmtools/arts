# Options Editor for ARTS Enum Option Groups

## Overview
Added a specialized editor for ARTS option group enums that provides a dropdown menu for selecting from available enum values.

## Implementation

### Files Created/Modified

1. **`python/src/pyarts3/gui/edit/Options.py`** (NEW)
   - Dropdown selector for enum option values
   - Displays type name and current value
   - Scrollable combo box for long option lists
   - Type-preserving selection

2. **`python/src/pyarts3/gui/edit/__init__.py`** (MODIFIED)
   - Added Options module import
   - Added dispatcher routing for option groups:
     ```python
     option_groups = arts.globals.option_groups()
     if type_name in option_groups:
         return Options.edit(value, parent=parent)
     ```

## How It Works

### Option Groups
ARTS defines 46 option group enum types accessible via:
```python
import pyarts3.arts as arts
option_groups = arts.globals.option_groups()
# Returns list like:
# ['ray_path_observer_agendaSetGeometricMaxStep',
#  'disort_settings_agenda_setup_layer_emission_type',
#  ...]
```

### Each Option Type
Each option group type (e.g., `ray_path_observer_agendaSetGeometricMaxStep`) provides:
- `get_options()` - Returns list of all valid enum values
- String representation via `str(option)` for display

### Editor Features
- **Type Detection**: Checks if `type(value).__name__` is in `option_groups`
- **Dropdown Menu**: QComboBox with all available options
- **Scrollable**: Set `maxVisibleItems=20` for long lists
- **Current Selection**: Pre-selects current value in dropdown
- **Type Preservation**: Returns the enum instance, not a string

### Example Types
```python
# ray_path_observer_agendaSetGeometricMaxStep
options = ['1/2', 'linear', '0']

# disort_settings_agenda_setup_layer_emission_type  
options = ['None', 'LinearInTau', 'LinearInTauNonLTE']
```

## Usage

### Programmatic
```python
from pyarts3.gui.edit import Options
import pyarts3.arts as arts

# Create an option instance
OptionType = arts.ray_path_observer_agendaSetGeometricMaxStep
current = OptionType.get_options()[1]  # 'linear'

# Edit it
new_value = Options.edit(current)
# Opens dialog with dropdown, returns selected option
```

### Through Dispatcher
```python
from pyarts3.gui.edit import edit
import pyarts3.arts as arts

option_value = arts.ray_path_observer_agendaSetGeometricMaxStep.get_options()[0]
result = edit.edit(option_value)
# Automatically routes to Options.edit()
```

### In Workspace/Results Tab
When a workspace variable or result is an option group enum, double-clicking automatically opens the Options editor.

## UI Layout

```
┌─────────────────────────────────────────┐
│ Edit ray_path_observer_agenda...       │
├─────────────────────────────────────────┤
│ Type: ray_path_observer_agendaSet...   │
│ Current value: linear                   │
│                                         │
│ Select option: [ linear       ▼ ]      │
│                 ├─ 1/2                  │
│                 ├─ linear  ◀── selected │
│                 └─ 0                    │
│                                         │
│               [ OK ] [ Cancel ]         │
└─────────────────────────────────────────┘
```

## Routing Priority

The dispatcher checks in this order:
1. Exact type module (e.g., `Vector.py` for `Vector`)
2. `ArrayOf*` → `ArrayOf.py`
3. **Option groups → `Options.py`** ✅ NEW
4. Has `__array__` → `NDarray.py`
5. Fallback → `Generic.py`

## Type Safety

✅ Input type: `ray_path_observer_agendaSetGeometricMaxStep`
✅ After edit: `ray_path_observer_agendaSetGeometricMaxStep`

The editor preserves the exact enum type by returning an instance from `get_options()`, not converting to string.

## Testing

Run `test_options_editor.py` to verify:
- ✅ 46 option groups detected
- ✅ Options retrieved via `get_options()`
- ✅ String representation works
- ✅ Type preserved after selection
- ✅ Dispatcher routes correctly
- ✅ Editor module imports successfully

Run `demo_options_editor.py` to see the actual UI.

## Benefits

1. **User-Friendly**: Dropdown prevents typos, shows all valid options
2. **Type-Safe**: Returns proper enum instance, not string
3. **Scalable**: Handles all 46 option groups automatically
4. **Discoverable**: Users can see all available options
5. **Consistent**: Same UI pattern as other editors
