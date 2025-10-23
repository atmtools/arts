# ARTS GUI Editors

This directory contains editor modules for ARTS workspace group types.

## Structure

Each editor is a separate Python module named after the workspace group type it edits:

```
edit/
├── __init__.py          # Auto-discovery and dispatch
├── Generic.py           # Fallback read-only viewer
├── Index.py             # Integer editor
├── Matrix.py            # 2D array editor
├── Numeric.py           # Float editor
├── Stokvec.py           # Stokes vector editor
├── String.py            # String editor
├── Vector.py            # 1D array editor
├── Vector2.py           # 2D vector editor
└── Vector3.py           # 3D vector editor
```

## Usage

### From Python code:

```python
import pyarts3.gui.edit as edit

# Edit any value - automatically dispatches to correct editor
vec = pyarts.arts.Vector([1, 2, 3])
new_vec = edit.edit(vec, parent=None)

# Or use specific editor directly
from pyarts3.gui.edit import Vector
new_vec = Vector.edit(vec)
```

### From GUI:

Double-click any item in the Simulation Settings, Results, or Additional Options lists to open its editor.

## Creating a New Editor

To add support for a new workspace group type:

### 1. Create the editor module

Create `edit/MyType.py`:

```python
"""Editor for MyType values."""

from PyQt5.QtWidgets import QDialog, QVBoxLayout, QDialogButtonBox
from PyQt5.QtCore import Qt

__all__ = ['edit']


def edit(value, parent=None):
    """
    Edit a MyType value.
    
    Parameters
    ----------
    value : MyType
        The value to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    MyType or None
        The edited value if accepted, None if cancelled
    """
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit MyType")
    
    layout = QVBoxLayout()
    
    # Add your custom editor widgets here
    # ...
    
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    if dialog.exec_() == QDialog.Accepted:
        # Return the edited value
        return modified_value
    return None
```

### 2. Import it in __init__.py

Add to `edit/__init__.py`:

```python
from . import MyType
```

That's it! The system will automatically discover and use your editor.

## Editor Guidelines

1. **Function signature**: `edit(value, parent=None)`
2. **Return value**: The edited value if OK clicked, `None` if cancelled
3. **Dialog**: Use `QDialog.exec_()` and check for `QDialog.Accepted`
4. **Buttons**: Include OK and Cancel buttons (use `QDialogButtonBox`)
5. **Type conversion**: Handle conversion to/from ARTS types if needed

## Fallback Behavior

If no specific editor exists for a type:
- The `Generic.py` editor is used
- Displays a read-only text view
- Always returns `None` (no editing)

## Examples

See existing editors for patterns:
- **Simple value**: `Numeric.py`, `Index.py`, `String.py`
- **Array editor**: `Vector.py`, `Matrix.py`
- **Component editor**: `Stokvec.py`, `Vector3.py`, `Vector2.py`
- **Read-only viewer**: `Generic.py`
