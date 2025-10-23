#!/usr/bin/env python3
"""
Visual demonstration of the Options editor.
This script shows what the editor looks like for option group enums.
"""

import sys
import pyarts3.arts as arts
from PyQt5.QtWidgets import QApplication

def test_options_editor():
    """Open the Options editor for a sample enum option."""
    
    # Create QApplication
    app = QApplication(sys.argv)
    
    # Get an option type
    OptionType = arts.ray_path_observer_agendaSetGeometricMaxStep
    options = OptionType.get_options()
    
    # Create an instance with a specific value
    current_value = options[1]  # 'linear'
    
    print("Opening Options editor...")
    print(f"Current value: {current_value}")
    print(f"Type: {type(current_value).__name__}")
    print(f"Available options: {[str(opt) for opt in options]}")
    print("\nThe editor will show:")
    print("  - Type name at the top")
    print("  - Current value")
    print("  - Dropdown menu with all options")
    print("  - OK/Cancel buttons")
    print("\nClose the dialog to continue...")
    
    # Open the editor
    from pyarts3.gui.edit import Options
    result = Options.edit(current_value)
    
    if result is not None:
        print(f"\n✓ User selected: {result}")
        print(f"  Type: {type(result).__name__}")
    else:
        print("\n✗ User cancelled")
    
    return result

if __name__ == "__main__":
    result = test_options_editor()
    sys.exit(0 if result is not None else 1)
