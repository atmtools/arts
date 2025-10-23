#!/usr/bin/env python3
"""
Example usage of workspace group editors in PlotGui.

This example demonstrates how double-clicking on workspace group values
opens appropriate editors for different types.
"""

import sys
import numpy as np

# Try to import pyarts3
try:
    from pyarts3.gui import start_gui
except ImportError:
    print("Error: Could not import pyarts3.gui")
    print("Make sure you're running from the ARTS build directory")
    sys.exit(1)


def example_plot_kwarg_func(**simulation_settings):
    """
    Example function that returns different workspace group types.
    
    This demonstrates various workspace groups that can be edited:
    - Numeric (float)
    - Index (int)
    - String (str)
    - Vector (1D numpy array)
    - Matrix (2D numpy array)
    """
    # Get simulation settings or use defaults
    frequency = simulation_settings.get('frequency', 1e9)
    temperature = simulation_settings.get('temperature', 300)
    
    # Return various workspace group types
    return {
        # Basic types
        'frequency': frequency,  # Numeric
        'n_levels': 50,  # Index
        'species': 'H2O-161',  # String
        
        # Array types
        'pressure_grid': np.logspace(5, -2, 50),  # Vector
        'temperature_field': np.random.randn(10, 10) + 273.15,  # Matrix
        
        # Fixed-size vectors
        'stokes_vector': np.array([1.0, 0.0, 0.0, 0.0]),  # Stokvec (4-element)
        'position': np.array([100e3, 200e3, 300e3]),  # Vector3 (3-element)
        
        # More complex examples
        'altitude': np.linspace(0, 100e3, 100),  # Vector
        'weights': np.ones(20) / 20,  # Vector
    }


def example_post_callback(fig, ax, **additional_options):
    """
    Post-processing callback to customize the plot.
    """
    ax.set_title(additional_options.get('title', 'ARTS Workspace Groups Demo'))
    ax.grid(True, alpha=0.3)
    if 'xlabel' in additional_options:
        ax.set_xlabel(additional_options['xlabel'])
    if 'ylabel' in additional_options:
        ax.set_ylabel(additional_options['ylabel'])


def example_plot_func(y, **kwargs):
    """
    Simple plotting function.
    """
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots()
    
    # Plot the primary data
    if isinstance(y, np.ndarray):
        if y.ndim == 1:
            ax.plot(y)
        elif y.ndim == 2:
            im = ax.imshow(y, aspect='auto', interpolation='nearest')
            fig.colorbar(im, ax=ax)
    
    return fig, ax


if __name__ == "__main__":
    # Initial simulation settings with various workspace groups
    initial_settings = {
        'frequency': 230e9,  # Numeric (Hz)
        'temperature': 280,   # Numeric (K)
        'zenith_angle': 45,   # Numeric (degrees)
    }
    
    print("=" * 60)
    print("ARTS Workspace Group Editors Demo")
    print("=" * 60)
    print()
    print("Double-click on any item in the three lists to edit:")
    print("  • Simulation Settings: Edit initial parameters")
    print("  • Results: Edit computed values")
    print("  • Additional Options: Edit plot options")
    print()
    print("Supported types:")
    print("  • Numeric (float): DoubleSpinBox editor")
    print("  • Index (int): SpinBox editor")
    print("  • String (str): LineEdit editor")
    print("  • Vector (1D array): Table editor")
    print("  • Matrix (2D array): Table editor")
    print("  • Stokvec/Vector3/Vector2: Component editors")
    print()
    print("Unsupported types show read-only view.")
    print("=" * 60)
    print()
    
    # Start the GUI
    start_gui(
        plot_kwarg_func=example_plot_kwarg_func,
        plot_func=example_plot_func,
        simulation_settings=initial_settings,
        post=example_post_callback
    )
