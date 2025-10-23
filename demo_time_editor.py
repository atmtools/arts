#!/usr/bin/env python3
"""Quick demo of Time editor features"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from PyQt5.QtWidgets import QApplication
from pyarts3 import arts
from pyarts3.gui.edit import edit
import datetime

app = QApplication(sys.argv)

# Demo: Edit current time
print("Opening Time editor with current datetime...")
print("Features:")
print("  - Calendar popup for date selection")
print("  - Arrow keys to increment/decrement values")
print("  - Format: yyyy-MM-dd HH:mm:ss.zzz (with milliseconds)")
print("  - Direct typing supported")
print()

t = arts.Time(datetime.datetime.now())
print(f"Current Time: {t.time}")

result = edit(t)

if result:
    print(f"Edited Time:  {result.time}")
    delta = result.time - t.time
    print(f"Time difference: {delta}")
else:
    print("Cancelled")
