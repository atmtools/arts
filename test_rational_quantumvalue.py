#!/usr/bin/env python3
"""Test Rational and QuantumValue editors"""

import sys
sys.path.insert(0, '/Users/richard/Work/arts/build/python/src')

from pyarts3 import arts
from pyarts3.gui.edit import edit

print("=" * 60)
print("Testing Rational editor...")
print("=" * 60)

r = arts.Rational(3, 2)
print(f"Original Rational:")
print(f"  String: {r}")
print(f"  n: {r.n}, d: {r.d}")

result_r = edit(r)

if result_r is not None:
    print(f"\nResult:")
    print(f"  String: {result_r}")
    print(f"  n: {result_r.n}, d: {result_r.d}")
    print(f"  Type preserved: {type(result_r) == type(r)}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("Testing QuantumValue editor (Rational type)...")
print("=" * 60)

qv_rational = arts.QuantumValue(arts.Rational(5, 2))
print(f"Original QuantumValue (Rational):")
print(f"  String: {qv_rational}")
print(f"  value: {qv_rational.value}")
print(f"  value type: {type(qv_rational.value).__name__}")

result_qv_r = edit(qv_rational)

if result_qv_r is not None:
    print(f"\nResult:")
    print(f"  String: {result_qv_r}")
    print(f"  value: {result_qv_r.value}")
    print(f"  value type: {type(result_qv_r.value).__name__}")
    print(f"  Type preserved: {type(result_qv_r) == type(qv_rational)}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("Testing QuantumValue editor (String type)...")
print("=" * 60)

qv_string = arts.QuantumValue('J')
print(f"Original QuantumValue (String):")
print(f"  String: {qv_string}")
print(f"  value: {qv_string.value}")
print(f"  value type: {type(qv_string.value).__name__}")

result_qv_s = edit(qv_string)

if result_qv_s is not None:
    print(f"\nResult:")
    print(f"  String: {result_qv_s}")
    print(f"  value: {result_qv_s.value}")
    print(f"  value type: {type(result_qv_s.value).__name__}")
    print(f"  Type preserved: {type(result_qv_s) == type(qv_string)}")
else:
    print("\nCancelled")

print("\n" + "=" * 60)
print("Testing full chain: QuantumIdentifier -> state -> QuantumUpperLower -> QuantumValue")
print("=" * 60)

qi = arts.QuantumIdentifier('O2-66 J 3 2')
print(f"Original QuantumIdentifier: {qi}")
print(f"  State J upper: {qi.state['J'].upper}")
print(f"  State J lower: {qi.state['J'].lower}")

result_qi = edit(qi)

if result_qi is not None:
    print(f"\nResult QuantumIdentifier: {result_qi}")
    if 'J' in result_qi.state:
        print(f"  State J upper: {result_qi.state['J'].upper}")
        print(f"  State J lower: {result_qi.state['J'].lower}")
else:
    print("\nCancelled")
