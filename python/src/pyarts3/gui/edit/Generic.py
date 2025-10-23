"""Generic editor for unsupported types with type-preserving eval option.

Behavior
--------
- If a specialized editor exists (a submodule in this package named after the
    value's type that exposes an `edit(value, parent=None)` function), delegate
    to that editor.
- Otherwise, fallback to a type-preserving eval UI where the user provides the
    constructor argument for the original type.
"""

import numpy as np
from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QTextEdit, QDialogButtonBox, QMessageBox
)
from PyQt5.QtCore import Qt
import sys

__all__ = ['edit']


def edit(value, parent=None):
    """
        Edit a value by preferring a specialized editor when available, with
        type-preserving eval as a fallback.

    The input box content is treated as an expression for the constructor of the
    original type. Concretely, this evaluates:

        eval(f"{TypeName}(<text>)")

    where TypeName is the original type's name bound in the eval environment.
    For example, for a String value, entering '"hello"' will call String("hello").

    If evaluation fails, an error is shown and the dialog stays open.

    Parameters
    ----------
    value : any
        The current value to view/edit.
    parent : QWidget, optional
        Parent widget for the dialog.

    Returns
    -------
    any or None
        The edited value if evaluation succeeds and OK is pressed; None if
        cancelled or evaluation fails.
    """
    original_type = type(value)
    type_name = original_type.__name__

    # Prefer a specialized editor if one exists in this package
    try:
        # Parent package module, e.g., 'pyarts3.gui.edit'
        parent_pkg_name = __name__.rsplit('.', 1)[0]
        parent_pkg = sys.modules.get(parent_pkg_name)
        if parent_pkg is None:
            parent_pkg = __import__(parent_pkg_name, fromlist=['*'])

        if hasattr(parent_pkg, type_name):
            submod = getattr(parent_pkg, type_name)
            # Avoid infinite recursion to ourselves and ensure edit() exists
            if getattr(submod, '__name__', None) != __name__ and hasattr(submod, 'edit'):
                return submod.edit(value, parent=parent)
    except Exception:
        # If any issue occurs during delegation discovery, silently fallback to eval UI
        pass

    dialog = QDialog(parent)
    dialog.setWindowTitle(f"Edit {type_name}")

    layout = QVBoxLayout()

    # Instruction
    info = QLabel(
        f"Type: {type_name}\n"
        "Enter a Python expression for the constructor argument.\n"
        f"Will evaluate: {type_name}(<expr>)\n"
        "Examples: 42, 3.14, 'text', [1, 2, 3], np.arange(5).tolist()"
    )
    info.setStyleSheet("color: #555;")
    layout.addWidget(info)

    # Editable text area initialized to repr(value)
    text_edit = QTextEdit()
    text_edit.setPlainText(repr(value))
    text_edit.setMinimumHeight(80)
    layout.addWidget(text_edit)

    # Buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    layout.addWidget(buttons)

    # Holder for the result
    result_holder = {"value": None}

    def on_accept():
        expr = text_edit.toPlainText().strip()
        # Build a safe-ish eval environment with the original type bound by name
        env = {
            '__builtins__': __builtins__,
            'np': np,
            'numpy': np,
            type_name: original_type,
        }
        try:
            constructed = eval(f"{type_name}({expr})", env, {})
            # Optional: ensure type preservation
            if not isinstance(constructed, original_type):
                # Try coercion as a fallback
                try:
                    constructed = original_type(constructed)
                except Exception:
                    raise TypeError(
                        f"Evaluation did not produce a {type_name}, got {type(constructed).__name__}")
            result_holder["value"] = constructed
            dialog.accept()
        except Exception as e:
            QMessageBox.warning(dialog, "Invalid edit",
                                f"Could not evaluate {type_name}(<expr>):\n{e}")
            # Keep dialog open

    buttons.accepted.connect(on_accept)
    buttons.rejected.connect(dialog.reject)

    dialog.setLayout(layout)

    if dialog.exec_() == QDialog.Accepted:
        return result_holder["value"]
    return None
