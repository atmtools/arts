"""
Editors for ARTS workspace group types.

This module provides editing capabilities for all workspace group types
using pyarts3.utils.builtin_groups() for type discovery.
"""

from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QTextEdit, 
                              QLineEdit, QSpinBox, QDoubleSpinBox, QPushButton,
                              QLabel, QDialogButtonBox, QFormLayout, QTableWidget,
                              QTableWidgetItem, QHeaderView, QTabWidget, QWidget,
                              QComboBox, QCheckBox, QScrollArea)
from PyQt5.QtCore import Qt
import numpy as np
from pyarts3.utils import builtin_groups


class BaseEditor(QDialog):
    """Base class for workspace group editors."""
    
    def __init__(self, value, parent=None):
        super().__init__(parent)
        self.original_value = value
        self.modified_value = value
        self.init_ui()
    
    def init_ui(self):
        """Initialize the UI. Override in subclasses."""
        layout = QVBoxLayout()
        
        # Add content area - override in subclass
        self.setup_content(layout)
        
        # Add buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        
        self.setLayout(layout)
    
    def setup_content(self, layout):
        """Setup the content area. Override in subclasses."""
        pass
    
    def get_value(self):
        """Get the modified value. Override in subclasses."""
        return self.modified_value


class EditGeneric(BaseEditor):
    """Generic text-based editor for unknown types."""
    
    def __init__(self, value, parent=None):
        super().__init__(value, parent)
        self.setWindowTitle(f"View {type(value).__name__}")
    
    def setup_content(self, layout):
        layout.addWidget(QLabel(f"Type: {type(self.original_value).__name__}"))
        layout.addWidget(QLabel("(Read-only view)"))
        
        self.text_edit = QTextEdit()
        self.text_edit.setPlainText(str(self.original_value))
        self.text_edit.setReadOnly(True)
        layout.addWidget(self.text_edit)


class EditNumeric(BaseEditor):
    """Editor for Numeric (floating point) values."""
    
    def __init__(self, value, parent=None):
        super().__init__(value, parent)
        self.setWindowTitle("Edit Numeric")
    
    def setup_content(self, layout):
        form = QFormLayout()
        
        self.spin_box = QDoubleSpinBox()
        self.spin_box.setRange(-1e308, 1e308)
        self.spin_box.setDecimals(10)
        self.spin_box.setValue(float(self.original_value))
        form.addRow("Value:", self.spin_box)
        
        layout.addLayout(form)
    
    def get_value(self):
        return self.spin_box.value()


class EditIndex(BaseEditor):
    """Editor for Index (integer) values."""
    
    def __init__(self, value, parent=None):
        super().__init__(value, parent)
        self.setWindowTitle("Edit Index")
    
    def setup_content(self, layout):
        form = QFormLayout()
        
        self.spin_box = QSpinBox()
        self.spin_box.setRange(-2147483648, 2147483647)
        self.spin_box.setValue(int(self.original_value))
        form.addRow("Value:", self.spin_box)
        
        layout.addLayout(form)
    
    def get_value(self):
        return self.spin_box.value()


class EditString(BaseEditor):
    """Editor for String values."""
    
    def __init__(self, value, parent=None):
        super().__init__(value, parent)
        self.setWindowTitle("Edit String")
    
    def setup_content(self, layout):
        form = QFormLayout()
        
        self.line_edit = QLineEdit()
        self.line_edit.setText(str(self.original_value))
        form.addRow("Value:", self.line_edit)
        
        layout.addLayout(form)
    
    def get_value(self):
        return self.line_edit.text()


class EditVector(BaseEditor):
    """Editor for Vector (1D array) values."""
    
    def __init__(self, value, parent=None):
        super().__init__(value, parent)
        self.setWindowTitle("Edit Vector")
        self.resize(600, 400)
    
    def setup_content(self, layout):
        # Convert to numpy array if needed
        if hasattr(self.original_value, '__array__'):
            self.data = np.array(self.original_value)
        else:
            self.data = np.array(self.original_value)
        
        # Info label
        info = QLabel(f"Size: {len(self.data)}")
        layout.addWidget(info)
        
        # Table
        self.table = QTableWidget()
        self.table.setColumnCount(2)
        self.table.setHorizontalHeaderLabels(["Index", "Value"])
        self.table.horizontalHeader().setStretchLastSection(True)
        
        self.table.setRowCount(len(self.data))
        for i, val in enumerate(self.data):
            # Index (read-only)
            idx_item = QTableWidgetItem(str(i))
            idx_item.setFlags(idx_item.flags() & ~Qt.ItemIsEditable)
            self.table.setItem(i, 0, idx_item)
            
            # Value (editable)
            val_item = QTableWidgetItem(f"{val:.10g}")
            self.table.setItem(i, 1, val_item)
        
        layout.addWidget(self.table)
    
    def get_value(self):
        result = np.zeros(len(self.data))
        for i in range(self.table.rowCount()):
            try:
                result[i] = float(self.table.item(i, 1).text())
            except (ValueError, AttributeError):
                result[i] = self.data[i]  # Keep original on error
        return result


class EditMatrix(BaseEditor):
    """Editor for Matrix (2D array) values."""
    
    def __init__(self, value, parent=None):
        super().__init__(value, parent)
        self.setWindowTitle("Edit Matrix")
        self.resize(800, 600)
    
    def setup_content(self, layout):
        # Convert to numpy array if needed
        if hasattr(self.original_value, '__array__'):
            self.data = np.array(self.original_value)
        else:
            self.data = np.array(self.original_value)
        
        # Info label
        info = QLabel(f"Shape: {self.data.shape[0]} × {self.data.shape[1]}")
        layout.addWidget(info)
        
        # Table
        self.table = QTableWidget()
        self.table.setRowCount(self.data.shape[0])
        self.table.setColumnCount(self.data.shape[1])
        
        # Set headers
        self.table.setHorizontalHeaderLabels([str(i) for i in range(self.data.shape[1])])
        self.table.setVerticalHeaderLabels([str(i) for i in range(self.data.shape[0])])
        
        # Fill data
        for i in range(self.data.shape[0]):
            for j in range(self.data.shape[1]):
                item = QTableWidgetItem(f"{self.data[i, j]:.6g}")
                self.table.setItem(i, j, item)
        
        # Make table scrollable
        scroll = QScrollArea()
        scroll.setWidget(self.table)
        scroll.setWidgetResizable(True)
        layout.addWidget(scroll)
    
    def get_value(self):
        result = np.zeros(self.data.shape)
        for i in range(self.data.shape[0]):
            for j in range(self.data.shape[1]):
                try:
                    result[i, j] = float(self.table.item(i, j).text())
                except (ValueError, AttributeError):
                    result[i, j] = self.data[i, j]  # Keep original on error
        return result


class EditStokvec(BaseEditor):
    """Editor for fixed-size vectors (Stokvec, Vector3, Vector2)."""
    
    def __init__(self, value, parent=None):
        super().__init__(value, parent)
        self.setWindowTitle("Edit Vector")
    
    def setup_content(self, layout):
        # Convert to numpy array if needed
        if hasattr(self.original_value, '__array__'):
            self.data = np.array(self.original_value)
        else:
            self.data = np.array(self.original_value)
        
        form = QFormLayout()
        self.spin_boxes = []
        
        labels = {
            4: ["I", "Q", "U", "V"],  # Stokvec
            3: ["X", "Y", "Z"],        # Vector3
            2: ["X", "Y"],              # Vector2
        }
        
        component_labels = labels.get(len(self.data), [str(i) for i in range(len(self.data))])
        
        for i, (label, val) in enumerate(zip(component_labels, self.data)):
            spin = QDoubleSpinBox()
            spin.setRange(-1e308, 1e308)
            spin.setDecimals(10)
            spin.setValue(float(val))
            form.addRow(f"{label}:", spin)
            self.spin_boxes.append(spin)
        
        layout.addLayout(form)
    
    def get_value(self):
        return np.array([spin.value() for spin in self.spin_boxes])


# Alias for Vector3 and Vector2 using the same implementation
EditVector3 = EditStokvec
EditVector2 = EditStokvec


def open_workspace_group_editor(value, parent=None):
    """
    Open an appropriate editor for the given workspace group value.
    
    Uses the naming convention Edit{TypeName} to find the editor class.
    Falls back to EditGeneric for types without specific editors.
    
    Parameters
    ----------
    value : any
        The value to edit
    parent : QWidget, optional
        Parent widget
    
    Returns
    -------
    modified_value or None
        The modified value if accepted, None if cancelled
    """
    type_name = type(value).__name__
    editor_name = f"Edit{type_name}"
    
    # Try to get the editor class from this module's globals
    editor_class = globals().get(editor_name)
    
    if editor_class is None:
        # Fall back to generic editor
        editor_class = EditGeneric
    
    editor = editor_class(value, parent)
    if editor.exec_() == QDialog.Accepted:
        return editor.get_value()
    return None
