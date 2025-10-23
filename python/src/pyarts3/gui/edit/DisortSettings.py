"""Editor for DisortSettings type."""

from PyQt5.QtWidgets import (
    QDialog, QVBoxLayout, QTableWidget, QTableWidgetItem,
    QDialogButtonBox, QHeaderView, QLabel, QScrollArea
)
from PyQt5.QtCore import Qt


def edit(value, parent=None):
    """
    Edit a DisortSettings object.
    
    DisortSettings contains configuration for DISORT radiative transfer:
    - altitude_grid: DescendingGrid
    - frequency_grid: AscendingGrid
    - quadrature_dimension: int
    - fourier_mode_dimension: int
    - legendre_polynomial_dimension: int
    - optical_thicknesses: Matrix
    - single_scattering_albedo: Matrix
    - fractional_scattering: Matrix
    - legendre_coefficients: Tensor3
    - solar_source: Vector
    - solar_zenith_angle: Vector
    - solar_azimuth_angle: Vector
    - positive_boundary_condition: Tensor3
    - negative_boundary_condition: Tensor3
    - source_polynomial: Tensor3
    - bidirectional_reflectance_distribution_functions: MatrixOfDisortBDRF
    
    Parameters
    ----------
    value : DisortSettings
        The DisortSettings object to edit
    parent : QWidget, optional
        Parent widget
        
    Returns
    -------
    DisortSettings or None
        Modified DisortSettings if OK clicked, None if cancelled
    """
    from pyarts3 import arts
    from pyarts3.gui.edit import edit as dispatch_edit
    
    dialog = QDialog(parent)
    dialog.setWindowTitle("Edit DisortSettings")
    dialog.setMinimumWidth(800)
    dialog.setMinimumHeight(600)
    
    layout = QVBoxLayout(dialog)
    
    layout.addWidget(QLabel("DISORT radiative transfer settings. Double-click Value to edit."))
    
    # Use scrollable area for many fields
    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    
    table = QTableWidget()
    table.setColumnCount(3)
    table.setHorizontalHeaderLabels(["Field", "Type", "Value"])
    table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
    
    fields = [
        ("altitude_grid", value.altitude_grid, "DescendingGrid"),
        ("frequency_grid", value.frequency_grid, "AscendingGrid"),
        ("quadrature_dimension", value.quadrature_dimension, "int"),
        ("fourier_mode_dimension", value.fourier_mode_dimension, "int"),
        ("legendre_polynomial_dimension", value.legendre_polynomial_dimension, "int"),
        ("optical_thicknesses", value.optical_thicknesses, "Matrix"),
        ("single_scattering_albedo", value.single_scattering_albedo, "Matrix"),
        ("fractional_scattering", value.fractional_scattering, "Matrix"),
        ("legendre_coefficients", value.legendre_coefficients, "Tensor3"),
        ("solar_source", value.solar_source, "Vector"),
        ("solar_zenith_angle", value.solar_zenith_angle, "Vector"),
        ("solar_azimuth_angle", value.solar_azimuth_angle, "Vector"),
        ("positive_boundary_condition", value.positive_boundary_condition, "Tensor3"),
        ("negative_boundary_condition", value.negative_boundary_condition, "Tensor3"),
        ("source_polynomial", value.source_polynomial, "Tensor3"),
        ("bidirectional_reflectance_distribution_functions", 
         value.bidirectional_reflectance_distribution_functions, "MatrixOfDisortBDRF"),
    ]
    
    table.setRowCount(len(fields))
    current_values = {}
    
    def refresh_table():
        for row, (field_name, field_value, field_type) in enumerate(fields):
            # Field name
            item = QTableWidgetItem(field_name)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 0, item)
            
            # Type
            item = QTableWidgetItem(field_type)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 1, item)
            
            # Value
            if field_name in current_values:
                val = current_values[field_name]
            else:
                val = field_value
                current_values[field_name] = val
            
            value_str = str(val)
            if len(value_str) > 100:
                value_str = value_str[:100] + "..."
            item = QTableWidgetItem(value_str)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            table.setItem(row, 2, item)
    
    refresh_table()
    
    def on_cell_double_clicked(row, col):
        if col != 2:  # Only edit Value column
            return
        
        field_name = fields[row][0]
        current_val = current_values[field_name]
        
        result = dispatch_edit(current_val, parent=dialog)
        if result is not None:
            current_values[field_name] = result
            refresh_table()
    
    table.cellDoubleClicked.connect(on_cell_double_clicked)
    scroll.setWidget(table)
    layout.addWidget(scroll)
    
    # OK/Cancel buttons
    buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    
    if dialog.exec_() == QDialog.Accepted:
        result = arts.DisortSettings()
        for field_name, _, _ in fields:
            setattr(result, field_name, current_values[field_name])
        return result
    
    return None
