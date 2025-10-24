Editor modules for GUI integration
==================================

Overview
--------

The ARTS Python GUI provides per-type editor modules for workspace
group values. These modules live in ``python/src/pyarts3/gui/edit/`` and
are used by the GUI (and can also be called directly) to edit values
in a small, focused dialog.

The goal is to give developers a clear, consistent pattern for adding or
adjusting editors as the set of supported ARTS types grows.

Design Pattern
--------------

Each editor module is named after the ARTS workspace type it handles
and exports a single function ``edit(value, parent=None)``:

.. code-block:: python

    __all__ = ["edit"]

    def edit(value, parent=None):
        """Edit a <WorkspaceGroup> value.

        Parameters
        ----------
        value : ~pyarts3.arts.<WorkspaceGroup>
            The value to edit
        parent : QWidget, optional
            Parent widget for the dialog

        Returns
        -------
        <WorkspaceGroup> or None
            The edited value if accepted, None if cancelled
        """
        # Build QDialog UI and return new value on accept

Dispatch and Fallback
---------------------

- The dispatcher ``pyarts3.gui.edit.edit(value, parent)`` looks up an
  editor submodule matching ``type(value).__name__`` and calls its
  ``edit`` function if available.
- If no specialized editor is found, the ``Generic`` editor is used as a
  fallback. It offers a type-preserving eval UI where the user enters an
  expression that is passed to the original type's constructor
  (``TypeName(<expr>)``). This makes it possible to edit uncommon types without
  writing a dedicated editor.
- For safety and consistency, the ``Generic`` editor first attempts to delegate
  to a specialized editor module if one exists; the eval UI is strictly a
  fallback.

Editor Guidelines
-----------------

- Signature: ``edit(value, parent=None)`` returning the edited value
  or ``None``.
- Dialogs should use ``QDialog.exec_()`` and return values only on
  ``QDialog.Accepted``.
- Preserve ARTS types. If the editor creates intermediate Python lists or numpy
  arrays, coerce back to the original type, e.g., ``type(value)(data)``.
- Keep the UI minimal and focused on the specific type.
- Provide friendly error messages and avoid crashing on invalid input.
- Respect the optional ``parent`` QWidget for modality/ownership.

Adding a new editor
-------------------

1. Create editor module at ``pyarts3/gui/edit/MyType.py``.
2. Implement a small ``QDialog`` to edit the value.
3. Reconstruct the type on save, e.g., ``type(value)(payload)``.
4. Add ``from . import MyType`` to ``pyarts3/gui/edit/__init__.py``:

   .. code-block:: python

       from . import MyType

5. Verify via GUI double-click or by calling the editor in Python:

   .. code-block:: python

       import pyarts3.gui.edit as edit
       new_value = edit.edit(old_value)

Examples and References
-----------------------

- Simple scalars: ``Numeric.py``, ``Index.py``, ``String.py``
- Arrays and matrices: ``Vector.py``, ``Matrix.py``
- Fixed-size vectors: ``Vector2.py``, ``Vector3.py``, ``Stokvec.py``
- Fallback: ``Generic.py`` (type-preserving eval)

These examples illustrate type preservation, validation, and concise UI design.
