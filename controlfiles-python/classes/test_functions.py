#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 09:59:42 2022

@author: u237023
"""

import os

def io(x, fname="tmp.xml", delete=False):
    try:
        x.savexml(fname)
        x.readxml(fname)
    finally:
        if os.path.exists(fname) and delete:
            os.remove(fname)

def array(x):
    assert len(x) == 1
    x.append(x[0])
    assert len(x) == 2
    x[0] = x[1]
    x.pop()
    assert len(x) == 1
