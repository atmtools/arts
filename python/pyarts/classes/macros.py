#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 08:40:53 2021

@author: u237023
"""

import ctypes as c
from pyarts.workspace.api import arts_api as lib


def BasicInterfaceCAPI(lib, type):
    exec(f"""
lib.create{type}.restype = c.c_void_p
lib.create{type}.argtypes = []

lib.delete{type}.restype = None
lib.delete{type}.argtypes = [c.c_void_p]

lib.print{type}.restype = None
lib.print{type}.argtypes = [c.c_void_p]
""")


def EnumMacroInterfaceCAP(lib, type):
        exec(f"""
lib.get{type}.restype = c.c_void_p
lib.get{type}.argtypes = [c.c_void_p]

lib.delete{type}.restype = c.c_int
lib.delete{type}.argtypes = [c.c_void_p, c.c_char_p]
""")

def VoidStructGetterCAPI(lib, type, value):
    exec(f"""
lib.get{value}{type}.restype = c.c_void_p
lib.get{value}{type}.argtypes = [c.c_void_p]
""")
