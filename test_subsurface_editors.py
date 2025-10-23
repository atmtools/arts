#!/usr/bin/env python3
"""Test SubsurfacePoint, SubsurfaceField, and SubsurfaceData editors"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from PyQt5.QtWidgets import QApplication
from pyarts3 import arts
from pyarts3.gui.edit import edit

# Create QApplication
app = QApplication(sys.argv)

print("=" * 60)
print("Testing SubsurfacePoint editor...")
print("=" * 60)

ssp = arts.SubsurfacePoint()
ssp.temperature = 280.0
ssp.density = 1000.0  # kg/m^3

print(f"Original SubsurfacePoint:")
print(f"  temperature: {ssp.temperature} K")
print(f"  density: {ssp.density} kg/m^3")
print(f"  keys: {list(ssp.keys())}")

result_ssp = edit(ssp)

if result_ssp is not None:
    print(f"\nResult SubsurfacePoint:")
    print(f"  temperature: {result_ssp.temperature} K")
    print(f"  density: {result_ssp.density} kg/m^3")
    print(f"  keys: {list(result_ssp.keys())}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("Testing SubsurfaceField editor...")
print("=" * 60)

ssf = arts.SubsurfaceField(5000.0)  # 5000 m bottom depth

print(f"Original SubsurfaceField:")
print(f"  bottom_depth: {ssf.bottom_depth} m")
print(f"  keys: {list(ssf.keys())}")

result_ssf = edit(ssf)

if result_ssf is not None:
    print(f"\nResult SubsurfaceField:")
    print(f"  bottom_depth: {result_ssf.bottom_depth} m")
    print(f"  keys: {list(result_ssf.keys())}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("Testing SubsurfaceData editor...")
print("=" * 60)

ssd = arts.SubsurfaceData(arts.Numeric(300.0))
print(f"Original SubsurfaceData:")
print(f"  data_type: {ssd.data_type}")
print(f"  alt_low: {ssd.alt_low}")
print(f"  alt_upp: {ssd.alt_upp}")
print(f"  lat_low: {ssd.lat_low}")
print(f"  lat_upp: {ssd.lat_upp}")
print(f"  lon_low: {ssd.lon_low}")
print(f"  lon_upp: {ssd.lon_upp}")

result_ssd = edit(ssd)

if result_ssd is not None:
    print(f"\nResult SubsurfaceData:")
    print(f"  data_type: {result_ssd.data_type}")
    print(f"  alt_low: {result_ssd.alt_low}")
    print(f"  alt_upp: {result_ssd.alt_upp}")
    print(f"  lat_low: {result_ssd.lat_low}")
    print(f"  lat_upp: {result_ssd.lat_upp}")
    print(f"  lon_low: {result_ssd.lon_low}")
    print(f"  lon_upp: {result_ssd.lon_upp}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("All tests completed")
print("=" * 60)
