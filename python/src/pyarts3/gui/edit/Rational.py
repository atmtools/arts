"""Editor for Rational type.

Rational has:
- n: numerator (int)
- d: denominator (int)
"""

from PyQt5.QtWidgets import (
    QApplication, QDialog, QVBoxLayout, QHBoxLayout, QLabel, QSpinBox,
    QDialogButtonBox
)

__all__ = ["edit"]


def edit(value, parent=None):
    """
    Edit a Rational value.
    
    Parameters
    ----------
    value : Rational
        The rational number to edit
    parent : QWidget, optional
        Parent widget for the dialog
    
    Returns
    -------
    Rational or None
        The edited value if accepted, None if cancelled
    """
    # Ensure a QApplication exists
    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit Rational")
    dialog.resize(400, 150)
    
    layout = QVBoxLayout()
    
    # Info
    info = QLabel(f"Current value: {value}")
    info.setStyleSheet("color: #555; font-weight: bold;")
    layout.addWidget(info)
    
    # Numerator
    n_layout = QHBoxLayout()
    n_layout.addWidget(QLabel("Numerator (n):"))
    n_spin = QSpinBox()
    n_spin.setMinimum(-1000000)
    n_spin.setMaximum(1000000)
    n_spin.setValue(value.n)
    n_layout.addWidget(n_spin)
    n_layout.addStretch()
    layout.addLayout(n_layout)
    
    # Denominator
    d_layout = QHBoxLayout()
    d_layout.addWidget(QLabel("Denominator (d):"))
    d_spin = QSpinBox()
    d_spin.setMinimum(-1000000)
    d_spin.setMaximum(1000000)
    d_spin.setValue(value.d)
    d_layout.addWidget(d_spin)
    d_layout.addStretch()
    layout.addLayout(d_layout)
    
    # Preview
    preview = QLabel(f"= {value.n}/{value.d}")
    preview.setStyleSheet("color: #666; font-size: 12px;")
    layout.addWidget(preview)
    
    def update_preview():
        n = n_spin.value()
        d = d_spin.value()
        preview.setText(f"= {n}/{d}")
    
    n_spin.valueChanged.connect(update_preview)
    d_spin.valueChanged.connect(update_preview)
    
    layout.addStretch()
    
    # Dialog buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    dialog.setLayout(layout)
    
    if dialog.exec_() == QDialog.Accepted:
        try:
            from pyarts3 import arts
            return arts.Rational(n_spin.value(), d_spin.value())
        except Exception as e:
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.warning(dialog, "Error", f"Cannot create Rational: {e}")
            return None
    return None
