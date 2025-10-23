#!/usr/bin/env python3
"""Test SurfacePoint, SurfaceField, and SurfaceData editors"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from PyQt5.QtWidgets import QApplication
from pyarts3 import arts
from pyarts3.gui.edit import edit

# Create QApplication
app = QApplication(sys.argv)

print("=" * 60)
print("Testing SurfacePoint editor...")
print("=" * 60)

sp = arts.SurfacePoint()
sp.temperature = 295.0
sp.elevation = 0.0
sp.normal = arts.Vector2([0.0, 1.0])

print(f"Original SurfacePoint:")
print(f"  temperature: {sp.temperature} K")
print(f"  elevation: {sp.elevation} m")
print(f"  normal: {sp.normal}")
print(f"  keys: {list(sp.keys())}")

result_sp = edit(sp)

if result_sp is not None:
    print(f"\nResult SurfacePoint:")
    print(f"  temperature: {result_sp.temperature} K")
    print(f"  elevation: {result_sp.elevation} m")
    print(f"  normal: {result_sp.normal}")
    print(f"  keys: {list(result_sp.keys())}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("Testing SurfaceField editor...")
print("=" * 60)

sf = arts.SurfaceField("Earth")

print(f"Original SurfaceField:")
print(f"  ellipsoid: {sf.ellipsoid}")
print(f"  keys: {list(sf.keys())}")

result_sf = edit(sf)

if result_sf is not None:
    print(f"\nResult SurfaceField:")
    print(f"  ellipsoid: {result_sf.ellipsoid}")
    print(f"  keys: {list(result_sf.keys())}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("Testing SurfaceData editor...")
print("=" * 60)

sd = arts.SurfaceData(arts.Numeric(295.0))
print(f"Original SurfaceData:")
print(f"  data_type: {sd.data_type}")
print(f"  lat_low: {sd.lat_low}")
print(f"  lat_upp: {sd.lat_upp}")
print(f"  lon_low: {sd.lon_low}")
print(f"  lon_upp: {sd.lon_upp}")

result_sd = edit(sd)

if result_sd is not None:
    print(f"\nResult SurfaceData:")
    print(f"  data_type: {result_sd.data_type}")
    print(f"  lat_low: {result_sd.lat_low}")
    print(f"  lat_upp: {result_sd.lat_upp}")
    print(f"  lon_low: {result_sd.lon_low}")
    print(f"  lon_upp: {result_sd.lon_upp}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("All tests completed")
print("=" * 60)
