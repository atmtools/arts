#!/usr/bin/env python3
"""Test SpeciesIsotope and QuantumUpperLower editors"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from pyarts3 import arts
from pyarts3.gui.edit import edit

print("=" * 60)
print("Testing SpeciesIsotope editor...")
print("=" * 60)

si = arts.SpeciesIsotope('O2-66')
print(f"Original SpeciesIsotope:")
print(f"  String: {si}")
print(f"  spec: {si.spec}")
print(f"  isotname: {si.isotname}")
print(f"  mass: {si.mass}")
print(f"  gi: {si.gi}")

result_si = edit(si)

if result_si is not None:
    print(f"\nResult:")
    print(f"  String: {result_si}")
    print(f"  spec: {result_si.spec}")
    print(f"  isotname: {result_si.isotname}")
    print(f"  Type preserved: {type(result_si) == type(si)}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("Testing QuantumUpperLower editor...")
print("=" * 60)

qi = arts.QuantumIdentifier('O2-66 J 1 2')
qul = qi.state['J']

print(f"Original QuantumUpperLower:")
print(f"  String: {qul}")
print(f"  upper: {qul.upper}")
print(f"  lower: {qul.lower}")

result_qul = edit(qul)

if result_qul is not None:
    print(f"\nResult:")
    print(f"  String: {result_qul}")
    print(f"  upper: {result_qul.upper}")
    print(f"  lower: {result_qul.lower}")
    print(f"  Type preserved: {type(result_qul) == type(qul)}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("Testing nested editing (QuantumIdentifier -> QuantumState -> QuantumUpperLower)")
print("=" * 60)

qi2 = arts.QuantumIdentifier('H2O J 1 2 Ka 3 4')
print(f"Original QuantumIdentifier: {qi2}")
print(f"  State: {qi2.state}")

result_qi = edit(qi2)

if result_qi is not None:
    print(f"\nResult QuantumIdentifier: {result_qi}")
    print(f"  State: {result_qi.state}")
else:
    print("\nCancelled")
