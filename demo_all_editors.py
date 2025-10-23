#!/usr/bin/env python3
"""Launch any of the 8 new editor demos interactively."""

import sys
import os

demos = {
    '1': ('Sun', 'demo_sun_editor.py', 'Solar properties editor'),
    '2': ('PropagationPathPoint', 'demo_propagation_path_point_editor.py', 'Ray path point editor'),
    '3': ('SpeciesTag', 'demo_species_tag_editor.py', 'Species identification editor'),
    '4': ('SensorObsel', 'demo_sensor_obsel_editor.py', 'Sensor observation element editor'),
    '5': ('JacobianTargets', 'demo_jacobian_targets_editor.py', 'Jacobian targets editor'),
    '6': ('DisortSettings', 'demo_disort_settings_editor.py', 'DISORT RT settings editor (16 fields!)'),
    '7': ('AbsorptionBand', 'demo_absorption_band_editor.py', 'Absorption band editor'),
    '8': ('AbsorptionLine', 'demo_absorption_line_editor.py', 'Spectroscopic line editor'),
}

def main():
    print("=" * 70)
    print("ARTS GUI EDITOR DEMOS - 8 New Priority Types")
    print("=" * 70)
    print("\nSelect a demo to run:\n")
    
    for key, (name, script, desc) in demos.items():
        print(f"  {key}. {name:25s} - {desc}")
    
    print("\n  0. Exit")
    print("\n" + "=" * 70)
    
    while True:
        try:
            choice = input("\nEnter choice (0-8): ").strip()
            
            if choice == '0':
                print("Goodbye!")
                sys.exit(0)
            
            if choice in demos:
                name, script, desc = demos[choice]
                print(f"\nLaunching {name} demo...")
                print(f"Description: {desc}\n")
                os.system(f"python3 {script}")
                print(f"\n{name} demo closed.")
                
                # Ask if user wants to run another
                again = input("\nRun another demo? (y/n): ").strip().lower()
                if again != 'y':
                    print("Goodbye!")
                    sys.exit(0)
            else:
                print(f"Invalid choice: {choice}. Please enter 0-8.")
        
        except KeyboardInterrupt:
            print("\n\nGoodbye!")
            sys.exit(0)
        except Exception as e:
            print(f"Error: {e}")


if __name__ == '__main__':
    main()
