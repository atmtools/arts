#!/usr/bin/env python3
"""Test AtmPoint, AtmField, and AtmData editors"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from PyQt5.QtWidgets import QApplication
from pyarts3 import arts
from pyarts3.gui.edit import edit

# Create QApplication
app = QApplication(sys.argv)

print("=" * 60)
print("Testing AtmPoint editor...")
print("=" * 60)

ap = arts.AtmPoint()
ap.temperature = 250.0
ap.pressure = 101325.0
ap.mag = arts.Vector3([1.0, 2.0, 3.0])
ap.wind = arts.Vector3([10.0, 5.0, 0.0])

print(f"Original AtmPoint:")
print(f"  temperature: {ap.temperature} K")
print(f"  pressure: {ap.pressure} Pa")
print(f"  mag: {ap.mag}")
print(f"  wind: {ap.wind}")

result_ap = edit(ap)

if result_ap is not None:
    print(f"\nResult AtmPoint:")
    print(f"  temperature: {result_ap.temperature} K")
    print(f"  pressure: {result_ap.pressure} Pa")
    print(f"  mag: {result_ap.mag}")
    print(f"  wind: {result_ap.wind}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("Testing AtmField editor...")
print("=" * 60)

af = arts.AtmField()
af.top_of_atmosphere = 100000.0  # 100 km

print(f"Original AtmField:")
print(f"  top_of_atmosphere: {af.top_of_atmosphere} m")
print(f"  nlte: {len(af.nlte)} items")
print(f"  specs: {len(af.specs)} items")
print(f"  isots: {len(af.isots)} items")
print(f"  other: {len(af.other)} items")
print(f"  ssprops: {len(af.ssprops)} items")

result_af = edit(af)

if result_af is not None:
    print(f"\nResult AtmField:")
    print(f"  top_of_atmosphere: {result_af.top_of_atmosphere} m")
    print(f"  nlte: {len(result_af.nlte)} items")
    print(f"  specs: {len(result_af.specs)} items")
    print(f"  isots: {len(result_af.isots)} items")
    print(f"  other: {len(result_af.other)} items")
    print(f"  ssprops: {len(result_af.ssprops)} items")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("Testing AtmData editor...")
print("=" * 60)

ad = arts.AtmData()
print(f"Original AtmData:")
print(f"  data_type: {ad.data_type}")
print(f"  alt_low: {ad.alt_low}")
print(f"  alt_upp: {ad.alt_upp}")
print(f"  lat_low: {ad.lat_low}")
print(f"  lat_upp: {ad.lat_upp}")
print(f"  lon_low: {ad.lon_low}")
print(f"  lon_upp: {ad.lon_upp}")

result_ad = edit(ad)

if result_ad is not None:
    print(f"\nResult AtmData:")
    print(f"  data_type: {result_ad.data_type}")
    print(f"  alt_low: {result_ad.alt_low}")
    print(f"  alt_upp: {result_ad.alt_upp}")
    print(f"  lat_low: {result_ad.lat_low}")
    print(f"  lat_upp: {result_ad.lat_upp}")
    print(f"  lon_low: {result_ad.lon_low}")
    print(f"  lon_upp: {result_ad.lon_upp}")
else:
    print("\nCancelled")
