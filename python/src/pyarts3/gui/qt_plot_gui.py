"""
Minimal Qt GUI starter for ARTS plots integration.

Usage:
    from pyarts3.qt_plot_gui import start_plot_gui
    start_plot_gui(plot_kwarg_func)

Where plot_kwarg_func is a callable (with __call__(self, **kwargs)) that returns a dict of kwargs for plots.plot().
"""
import os
import sys
import platform as _platform

# macOS-specific: ensure platform.mac_ver() reports a version for Qt/Matplotlib compatibility
if sys.platform == "darwin":
    os.environ.setdefault("SYSTEM_VERSION_COMPAT", "1")
    os.environ.setdefault("MPLBACKEND", "Qt5Agg")
    os.environ.setdefault("QT_API", "PyQt5")
    # Work around matplotlib parsing of platform.mac_ver() returning '' on some setups
    try:
        if not _platform.mac_ver()[0]:
            _platform.mac_ver = lambda: ("10.16", ("", "", ""), "")  # type: ignore
    except Exception:
        pass

import matplotlib
matplotlib.use("Qt5Agg", force=True)
 
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from PyQt5.QtWidgets import (QApplication, QVBoxLayout, QHBoxLayout, QPushButton, 
                              QWidget, QListWidget, QLabel, QDialog, QLineEdit, 
                              QDialogButtonBox, QFormLayout, QListWidgetItem, QMenu, QAction,
                              QMessageBox, QTextEdit)
from PyQt5.QtCore import Qt
import numpy as np
from . import edit as editors


class PlotGui(QWidget):
    def __init__(self, plot_kwarg_func, plot_func, simulation_settings=None, post=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.plot_kwarg_func = plot_kwarg_func
        self.plot_func = plot_func
        self.simulation_settings = simulation_settings or {}
        self.results_kwargs = {}
        self.additional_options = {}
        self.last_added_option = None  # Track last added option for error recovery
        self.post = post
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("ARTS GUI")
        self.fig = Figure(figsize=(16, 12), constrained_layout=True)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setMinimumSize(800, 600)
        # Enable native matplotlib pan/zoom interactions
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self._panning = False
        self._pan_start = None
        self._zooming = False
        self._zoom_start = None
        self._zoom_rect = None
        self._zoom_axes = None
        self._original_limits = {}  # Store original axis limits for reset
        self.ax = None
        
        # Create main horizontal layout
        main_layout = QHBoxLayout()
        main_layout.setContentsMargins(10, 10, 10, 10)
        
        # Left side: canvas
        main_layout.addWidget(self.canvas, stretch=1)
        
        # Right side: control panel
        right_panel = QVBoxLayout()
        self.run_btn = QPushButton("Run")
        self.run_btn.clicked.connect(self.run_plot)
        right_panel.addWidget(self.run_btn)
        
        # Simulation Settings list
        sim_label = QLabel("Simulation Settings")
        right_panel.addWidget(sim_label)
        self.simulation_list = QListWidget()
        self.simulation_list.itemClicked.connect(lambda item: self.on_item_selected("simulation", item))
        self.simulation_list.itemDoubleClicked.connect(lambda item: self.on_item_double_clicked("simulation", item))
        right_panel.addWidget(self.simulation_list)
        
        # Results list
        results_label = QLabel("Results")
        right_panel.addWidget(results_label)
        self.results_list = QListWidget()
        self.results_list.itemClicked.connect(lambda item: self.on_item_selected("results", item))
        self.results_list.itemDoubleClicked.connect(lambda item: self.on_item_double_clicked("results", item))
        right_panel.addWidget(self.results_list)
        
        # Additional Options list
        options_label = QLabel("Additional Options")
        right_panel.addWidget(options_label)
        
        # List widget for options with context menu
        self.options_list = QListWidget()
        self.options_list.setMaximumHeight(200)
        self.options_list.setContextMenuPolicy(Qt.CustomContextMenu)
        self.options_list.customContextMenuRequested.connect(self.show_options_context_menu)
        self.options_list.itemDoubleClicked.connect(lambda item: self.on_item_double_clicked("options", item))
        self.options_list.itemClicked.connect(self.on_options_item_clicked)
        right_panel.addWidget(self.options_list)
        
        main_layout.addLayout(right_panel)
        self.setLayout(main_layout)

    def run_plot(self):
        try:
            # Clear figure completely
            self.fig.clear()
            self.ax = None
            self._original_limits.clear()

            # Populate Simulation Settings list
            self.simulation_list.clear()
            for key, value in self.simulation_settings.items():
                # Show only key and type name in list
                type_name = type(value).__name__
                self.simulation_list.addItem(f"{key}: <{type_name}>")

            # Step 1: Get Results kwargs from plot_kwarg_func, passing Simulation Settings
            self.results_kwargs = self.plot_kwarg_func(**self.simulation_settings)

            # Populate the Results list
            self.results_list.clear()
            for key, value in self.results_kwargs.items():
                # Show only key and type name in list
                type_name = type(value).__name__
                self.results_list.addItem(f"{key}: <{type_name}>")

            # Populate Additional Options list
            self.update_options_display()

            # Step 2: Merge Additional Options into Results (Additional Options override)
            plot_kwargs = {**self.results_kwargs, **self.additional_options}

            # Extract 'y' as first positional argument if it exists
            y_data = plot_kwargs.pop('y', None)
            # Ensure no 'post' is forwarded to the plotting function; GUI handles post-callback itself
            plot_kwargs.pop('post', None)

            # Step 3: Inject fig/ax
            plot_kwargs['fig'] = self.fig
            plot_kwargs['ax'] = None  # Let plot create axes as needed

            # Step 4: Call plot function with y as first arg (if exists), rest as kwargs
            if y_data is not None:
                _, self.ax = self.plot_func(y_data, **plot_kwargs)
            else:
                _, self.ax = self.plot_func(**plot_kwargs)

            # Step 5: Call post (if provided) inside the GUI, passing additional options
            if callable(self.post):
                self.post(self.fig, self.ax, **self.additional_options)

            # Remove any empty/unused axes that might have been created
            for ax in self.fig.get_axes():
                if not ax.has_data() and not ax.get_title() and not ax.get_xlabel() and not ax.get_ylabel():
                    self.fig.delaxes(ax)
            # Store original limits for all axes
            for ax in self.fig.get_axes():
                self._original_limits[ax] = {
                    'xlim': ax.get_xlim(),
                    'ylim': ax.get_ylim()
                }
            self.canvas.draw()

        except Exception as e:
            # Show error dialog with scrollable text
            error_msg = f"{type(e).__name__}: {str(e)}"
            
            # Extract last line for main display
            error_lines = error_msg.split('\n')
            last_line = error_lines[-1].strip() if error_lines else error_msg

            msg_box = QMessageBox(self)
            msg_box.setIcon(QMessageBox.Critical)
            msg_box.setWindowTitle("Plot Error")
            msg_box.setText(f"An error occurred while plotting:\n\n{last_line}")
            
            # Use detailed text which is automatically scrollable
            msg_box.setDetailedText(error_msg)
            msg_box.addButton(QMessageBox.Ok)
            
            # Set a reasonable size for the dialog
            msg_box.setStyleSheet("QTextEdit { min-width: 600px; min-height: 300px; }")
            msg_box.exec_()
    
    def on_item_double_clicked(self, list_type, item):
        """Handle double-click on an item from any of the three lists."""
        # Parse the item text to extract key and value
        text = item.text()
        
        # Check if it's the "+ Add option..." item
        if text.startswith("+ Add option"):
            self.add_option()
            return
        
        # Parse "key: value" or "key: value (type)" format
        if ": " not in text:
            return
        
        key, rest = text.split(": ", 1)
        
        # Get the actual value from the appropriate dict
        if list_type == "simulation":
            value = self.simulation_settings.get(key)
        elif list_type == "results":
            value = self.results_kwargs.get(key)
        elif list_type == "options":
            # For options list, the key is stored in UserRole
            key = item.data(Qt.UserRole)
            if key is None:  # "+ Add option..." item
                self.add_option()
                return
            value = self.additional_options.get(key)
        else:
            return
        
        if value is None:
            return
        
        # Open editor (will automatically detect type and fall back to generic viewer)
        new_value = editors.edit(value, self)
        
        # Update the value if changed
        if new_value is not None:
            if list_type == "simulation":
                self.simulation_settings[key] = new_value
                # Update display - show only type name
                self.simulation_list.clear()
                for k, v in self.simulation_settings.items():
                    type_name = type(v).__name__
                    self.simulation_list.addItem(f"{k}: <{type_name}>")
                # Re-run plot to regenerate results
                self.run_plot()
            elif list_type == "results":
                self.results_kwargs[key] = new_value
                # Update display - show only type name
                self.results_list.clear()
                for k, v in self.results_kwargs.items():
                    type_name = type(v).__name__
                    self.results_list.addItem(f"{k}: <{type_name}>")
                # Replot with updated results
                self.replot()
            elif list_type == "options":
                self.additional_options[key] = new_value
                self.update_options_display()
                # Replot with updated options
                self.replot()
    
    def show_options_context_menu(self, position):
        """Show context menu for Additional Options list."""
        menu = QMenu()
        
        # Always show "Add Option"
        add_action = QAction("Add Option...", self)
        add_action.triggered.connect(self.add_option)
        menu.addAction(add_action)
        
        # Show "Edit Option" and "Remove Option" only if a real item is selected
        item = self.options_list.itemAt(position)
        if item and item.data(Qt.UserRole) is not None:
            edit_action = QAction("Edit Option...", self)
            edit_action.triggered.connect(self.edit_option)
            menu.addAction(edit_action)
            
            remove_action = QAction("Remove Option", self)
            remove_action.triggered.connect(self.remove_option)
            menu.addAction(remove_action)
        
        # Show menu at cursor position
        menu.exec_(self.options_list.mapToGlobal(position))
    
    def update_options_display(self):
        """Refresh the Additional Options display."""
        self.options_list.clear()
        for key, value in self.additional_options.items():
            type_name = type(value).__name__
            item = QListWidgetItem(f"{key}: {value} ({type_name})")
            item.setData(Qt.UserRole, key)  # Store key for removal
            self.options_list.addItem(item)
        
        # Add a "+ Add option..." item at the end
        add_item = QListWidgetItem("+ Add option...")
        add_item.setData(Qt.UserRole, None)  # Mark as the add button
        add_item.setForeground(Qt.gray)
        self.options_list.addItem(add_item)
    
    def on_options_item_clicked(self, item):
        """Handle click on options list items."""
        # Check if it's the "+ Add option..." item
        if item.data(Qt.UserRole) is None:
            self.add_option()
    
    def add_option(self):
        """Show dialog to add a new Additional Option."""
        dialog = AddOptionDialog(self, existing_keys=self.additional_options.keys())
        if dialog.exec_() == QDialog.Accepted:
            key, value = dialog.get_values()
            if key:
                self.additional_options[key] = value
                self.last_added_option = key  # Track for error recovery
                self.update_options_display()
                # Replot with new options
                self.replot()
    
    def edit_option(self):
        """Edit an existing Additional Option."""
        current_item = self.options_list.currentItem()
        if current_item:
            key = current_item.data(Qt.UserRole)
            # Don't edit if it's the add button
            if key is not None and key in self.additional_options:
                current_value = self.additional_options[key]
                # Pass existing keys except the one being edited
                other_keys = [k for k in self.additional_options.keys() if k != key]
                dialog = AddOptionDialog(
                    self, 
                    existing_keys=other_keys, 
                    edit_mode=True, 
                    initial_key=key, 
                    initial_value=current_value
                )
                result = dialog.exec_()
                if dialog.deleted:
                    # Delete was clicked
                    del self.additional_options[key]
                    # Clear last_added_option if we deleted it
                    if self.last_added_option == key:
                        self.last_added_option = None
                    self.update_options_display()
                    # Replot after deletion
                    self.replot()
                elif result == QDialog.Accepted:
                    new_key, new_value = dialog.get_values()
                    if new_key:
                        # Remove old key if it changed
                        if new_key != key:
                            del self.additional_options[key]
                            # Update last_added_option if key changed
                            if self.last_added_option == key:
                                self.last_added_option = new_key
                        self.additional_options[new_key] = new_value
                        self.update_options_display()
                        # Replot with modified options
                        self.replot()
    
    def remove_option(self):
        """Remove the selected Additional Option."""
        current_item = self.options_list.currentItem()
        if current_item:
            key = current_item.data(Qt.UserRole)
            # Don't remove if it's the add button
            if key is not None and key in self.additional_options:
                del self.additional_options[key]
                # Clear last_added_option if we removed it
                if self.last_added_option == key:
                    self.last_added_option = None
                self.update_options_display()
                # Replot after removal
                self.replot()
    
    def replot(self):
        """Replot with current settings - similar to run_plot but without clearing simulation data."""
        if not self.results_kwargs:
            # No data to plot yet
            return
        
        try:
            # Clear figure
            self.fig.clear()
            self.ax = None
            self._original_limits.clear()
            
            # Step 2: Merge Additional Options into Results (Additional Options override)
            plot_kwargs = {**self.results_kwargs, **self.additional_options}

            # Extract 'y' as first positional argument if it exists
            y_data = plot_kwargs.pop('y', None)
            # Ensure no 'post' is forwarded to the plotting function; GUI handles post-callback itself
            plot_kwargs.pop('post', None)

            # Step 3: Inject fig/ax
            plot_kwargs['fig'] = self.fig
            plot_kwargs['ax'] = None  # Let plot create axes as needed

            # Step 4: Call plot function with y as first arg (if exists), rest as kwargs
            if y_data is not None:
                _, self.ax = self.plot_func(y_data, **plot_kwargs)
            else:
                _, self.ax = self.plot_func(**plot_kwargs)
            
            # Step 5: Call post (if provided) inside the GUI, passing additional options
            if callable(self.post):
                self.post(self.fig, self.ax, **self.additional_options)
            
            # Remove any empty/unused axes that might have been created
            for ax in self.fig.get_axes():
                if not ax.has_data() and not ax.get_title() and not ax.get_xlabel() and not ax.get_ylabel():
                    self.fig.delaxes(ax)
            # Store original limits for all axes
            for ax in self.fig.get_axes():
                self._original_limits[ax] = {
                    'xlim': ax.get_xlim(),
                    'ylim': ax.get_ylim()
                }
            self.canvas.draw()
            
        except Exception as e:
            # Show error dialog with scrollable text and option to delete or edit last addition
            error_msg = f"{type(e).__name__}: {str(e)}"
            
            # Extract last line for main display
            error_lines = error_msg.split('\n')
            last_line = error_lines[-1].strip() if error_lines else error_msg
            
            msg_box = QMessageBox(self)
            msg_box.setIcon(QMessageBox.Critical)
            msg_box.setWindowTitle("Plot Error")
            msg_box.setText(f"An error occurred while plotting:\n\n{last_line}")
            
            # Use detailed text which is automatically scrollable
            msg_box.setDetailedText(error_msg)
            
            # Only show "Delete" and "Edit" if we have a last added option
            if self.last_added_option and self.last_added_option in self.additional_options:
                msg_box.setInformativeText(f"Last addition: '{self.last_added_option}'")
                delete_btn = msg_box.addButton("Delete", QMessageBox.DestructiveRole)
                edit_btn = msg_box.addButton("Edit", QMessageBox.AcceptRole)
                msg_box.setDefaultButton(edit_btn)
            else:
                msg_box.addButton(QMessageBox.Ok)
            
            # Set a reasonable size for the dialog
            msg_box.setStyleSheet("QTextEdit { min-width: 600px; min-height: 300px; }")
            
            result = msg_box.exec_()
            
            # If user chose to delete or edit last addition
            if self.last_added_option and self.last_added_option in self.additional_options:
                clicked_button = msg_box.clickedButton()
                button_text = clicked_button.text()
                
                if button_text == "Delete":
                    # Delete the problematic option
                    del self.additional_options[self.last_added_option]
                    self.last_added_option = None
                    self.update_options_display()
                    # Try to replot without the problematic option
                    self.replot()
                elif button_text == "Edit":
                    # Open edit dialog for the problematic option
                    key = self.last_added_option
                    current_value = self.additional_options[key]
                    # Pass existing keys except the one being edited
                    other_keys = [k for k in self.additional_options.keys() if k != key]
                    dialog = AddOptionDialog(
                        self, 
                        existing_keys=other_keys, 
                        edit_mode=True, 
                        initial_key=key, 
                            initial_value=current_value,
                            allow_name_edit=True
                    )
                    edit_result = dialog.exec_()
                    if dialog.deleted:
                        # Delete was clicked in edit dialog
                        del self.additional_options[key]
                        self.last_added_option = None
                        self.update_options_display()
                        # Try to replot without the problematic option
                        self.replot()
                    elif edit_result == QDialog.Accepted:
                        new_key, new_value = dialog.get_values()
                        if new_key:
                            # Remove old key if it changed
                            if new_key != key:
                                del self.additional_options[key]
                                # Update last_added_option if key changed
                                self.last_added_option = new_key
                            self.additional_options[new_key] = new_value
                            self.update_options_display()
                            # Replot with modified options
                            self.replot()

    def on_scroll(self, event):
        """Handle mouse scroll for zooming."""
        if event.inaxes is None:
            return
        ax = event.inaxes
        # Get current limits
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        # Calculate zoom factor
        zoom_factor = 1.2 if event.button == 'down' else 0.8
        # Get mouse position in data coordinates
        xdata = event.xdata
        ydata = event.ydata
        # Calculate new limits centered on mouse position
        x_left = xdata - (xdata - xlim[0]) * zoom_factor
        x_right = xdata + (xlim[1] - xdata) * zoom_factor
        y_bottom = ydata - (ydata - ylim[0]) * zoom_factor
        y_top = ydata + (ylim[1] - ydata) * zoom_factor
        # Apply new limits
        ax.set_xlim([x_left, x_right])
        ax.set_ylim([y_bottom, y_top])
        self.canvas.draw_idle()

    def on_press(self, event):
        """Handle mouse button press for panning, zooming, and reset."""
        if event.inaxes is None:
            return
        # Double-click to restore original view
        if event.dblclick and event.inaxes in self._original_limits:
            ax = event.inaxes
            limits = self._original_limits[ax]
            ax.set_xlim(limits['xlim'])
            ax.set_ylim(limits['ylim'])
            self.canvas.draw_idle()
            return
        # Left-click for panning
        if event.button == 1:
            self._panning = True
            self._pan_start = (event.xdata, event.ydata)
            self._pan_axes = event.inaxes
        # Right-click for zoom-to-rectangle
        elif event.button == 3:
            self._zooming = True
            self._zoom_start = (event.xdata, event.ydata)
            self._zoom_axes = event.inaxes
            # Create a rectangle patch for visual feedback
            self._zoom_rect = Rectangle(
                (event.xdata, event.ydata), 0, 0,
                fill=False, edgecolor='red', linewidth=2, linestyle='--'
            )
            event.inaxes.add_patch(self._zoom_rect)
            self.canvas.draw_idle()

    def on_release(self, event):
        """Handle mouse button release."""
        # Handle zoom rectangle completion
        if self._zooming and self._zoom_rect is not None and event.inaxes == self._zoom_axes:
            # Remove the rectangle
            self._zoom_rect.remove()
            self._zoom_rect = None
            # Apply zoom if we have valid coordinates
            if (event.xdata is not None and event.ydata is not None and 
                self._zoom_start is not None):
                x0, y0 = self._zoom_start
                x1, y1 = event.xdata, event.ydata
                # Make sure we have a proper rectangle (not just a click)
                if abs(x1 - x0) > 0.01 * (self._zoom_axes.get_xlim()[1] - self._zoom_axes.get_xlim()[0]) and \
                   abs(y1 - y0) > 0.01 * (self._zoom_axes.get_ylim()[1] - self._zoom_axes.get_ylim()[0]):
                    # Set new limits (ensure proper ordering)
                    self._zoom_axes.set_xlim([min(x0, x1), max(x0, x1)])
                    self._zoom_axes.set_ylim([min(y0, y1), max(y0, y1)])
            self.canvas.draw_idle()
            self._zooming = False
            self._zoom_start = None
            self._zoom_axes = None
        # Handle panning release
        self._panning = False
        self._pan_start = None

    def on_motion(self, event):
        """Handle mouse motion for panning and zoom rectangle."""
        # Handle zoom rectangle drawing
        if self._zooming and self._zoom_rect is not None:
            if event.xdata is not None and event.ydata is not None and self._zoom_start is not None:
                x0, y0 = self._zoom_start
                width = event.xdata - x0
                height = event.ydata - y0
                self._zoom_rect.set_width(width)
                self._zoom_rect.set_height(height)
                self.canvas.draw_idle()
            return
        # Handle panning
        if not self._panning or event.inaxes != self._pan_axes:
            return
        if event.xdata is None or event.ydata is None or self._pan_start is None:
            return
        # Calculate delta
        dx = self._pan_start[0] - event.xdata
        dy = self._pan_start[1] - event.ydata
        # Update limits
        ax = self._pan_axes
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.set_xlim([xlim[0] + dx, xlim[1] + dx])
        ax.set_ylim([ylim[0] + dy, ylim[1] + dy])
        self.canvas.draw_idle()


class AddOptionDialog(QDialog):
    """Dialog for adding a new Additional Option with type inference."""
    
    def __init__(self, parent=None, existing_keys=None, edit_mode=False, initial_key=None, initial_value=None, allow_name_edit=False):
        super().__init__(parent)
        self.existing_keys = set(existing_keys) if existing_keys else set()
        self.edit_mode = edit_mode
        self.initial_key = initial_key
        self.allow_name_edit = allow_name_edit
        
        if edit_mode:
            self.setWindowTitle("Edit Additional Option")
        else:
            self.setWindowTitle("Add Additional Option")
        
        layout = QFormLayout()
        
        self.key_input = QLineEdit()
        self.key_input.setPlaceholderText("Option name")
        if edit_mode and initial_key:
            self.key_input.setText(initial_key)
            # Only gray out if not allowing name edit
            if not allow_name_edit:
                self.key_input.setEnabled(False)  # Gray out in edit mode
                self.key_input.setStyleSheet("QLineEdit { background-color: #f0f0f0; color: #888; }")
        self.key_input.textChanged.connect(self.check_duplicate)
        layout.addRow("Name:", self.key_input)
        
        self.warning_label = QLabel("")
        self.warning_label.setStyleSheet("QLabel { color: red; font-size: 10px; }")
        layout.addRow("", self.warning_label)
        
        self.value_input = QLineEdit()
        self.value_input.setPlaceholderText("Value (auto-typed; empty → None)")
        if edit_mode and initial_value is not None:
            self.value_input.setText(str(initial_value))
        layout.addRow("Value:", self.value_input)
        
        # Buttons
        if edit_mode:
            self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
            delete_btn = self.buttons.addButton("Delete", QDialogButtonBox.DestructiveRole)
            delete_btn.clicked.connect(self.delete_option)
        else:
            self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        
        self.ok_button = self.buttons.button(QDialogButtonBox.Ok)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addRow(self.buttons)
        
        self.setLayout(layout)
        self.deleted = False
        
        # Initial duplicate check if in edit mode
        if not edit_mode:
            self.check_duplicate()
    
    def check_duplicate(self):
        """Check if the key already exists and show warning."""
        key = self.key_input.text().strip()
        if key in self.existing_keys:
            self.warning_label.setText("Duplicate name")
            self.key_input.setStyleSheet("QLineEdit { border: 1px solid red; }")
            self.ok_button.setEnabled(False)
        else:
            self.warning_label.setText("")
            self.key_input.setStyleSheet("")
            self.ok_button.setEnabled(True)
    
    def delete_option(self):
        """Mark option for deletion and close dialog."""
        self.deleted = True
        self.reject()
    
    def get_values(self):
        """Return the key and inferred-type value."""
        if self.deleted:
            return None, None
        
        key = self.key_input.text().strip()
        value_str = self.value_input.text().strip()
        
        # Infer type
        value = self._infer_type(value_str)
        
        return key, value
    
    def _infer_type(self, value_str):
        """Infer type using eval with fallback to string; empty -> None.

        Behavior:
        - Empty string: return None
        - Non-empty: try eval(value_str); on any error, return the original string
        """
        # Empty string -> None
        if value_str == "":
            return None

        try:
            # Use eval to allow Python literals/expressions (e.g., lists, dicts, None, True, False)
            # Intentionally not restricting builtins here as this is user-entered config in a local GUI.
            return eval(value_str)
        except Exception:
            # Fallback to raw string if eval fails
            return value_str


def start_gui(plot_kwarg_func, plot_func, simulation_settings=None, post=None):
    app = QApplication.instance() or QApplication(sys.argv)
    gui = PlotGui(plot_kwarg_func, plot_func, simulation_settings=simulation_settings, post=post)
    gui.show()
    app.exec_()
