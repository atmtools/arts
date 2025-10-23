#!/usr/bin/env python3
"""Test Time editor"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from PyQt5.QtWidgets import QApplication
from pyarts3 import arts
from pyarts3.gui.edit import edit
import datetime

# Create QApplication
app = QApplication(sys.argv)

print("=" * 60)
print("Testing Time editor...")
print("=" * 60)

# Create a Time object with current datetime
now = datetime.datetime.now()
t = arts.Time(now)

print(f"Original Time:")
print(f"  Time object: {t}")
print(f"  .time attribute: {t.time}")
print(f"  Type: {type(t.time)}")

result_t = edit(t)

if result_t is not None:
    print(f"\nResult Time:")
    print(f"  Time object: {result_t}")
    print(f"  .time attribute: {result_t.time}")
    print(f"  Type: {type(result_t.time)}")
    print(f"\nTime changed: {t.time != result_t.time}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("Testing with specific datetime...")
print("=" * 60)

# Create with a specific datetime
specific_dt = datetime.datetime(2025, 1, 15, 14, 30, 45, 123456)
t2 = arts.Time(specific_dt)

print(f"Original Time:")
print(f"  Time object: {t2}")
print(f"  .time attribute: {t2.time}")

result_t2 = edit(t2)

if result_t2 is not None:
    print(f"\nResult Time:")
    print(f"  Time object: {result_t2}")
    print(f"  .time attribute: {result_t2.time}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("All tests completed")
print("=" * 60)
