#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 09:59:42 2022

@author: u237023
"""

import os
import uuid


def io(x, fname=None, delete=False):
    if fname is None:
        fname = str(uuid.uuid4()) + ".xml"
    try:
        x.savexml(fname)
        x.readxml(fname)
    finally:
        if os.path.exists(fname) and delete:
            os.remove(fname)
    return True

def array(x):
    assert len(x) == 1
    x.append(x[0])
    assert len(x) == 2
    x[0] = x[1]
    x.pop()
    assert len(x) == 1
    y = type(x)([x[0]])
    assert len(y) == 1
    return True

def array_of_array(x):
    assert len(x) > 0
    y = type(x)([[x[0][0]]])
    assert len(y) == 1
    assert len(y[0]) == 1
    return True

def shape_match(x, y):
    x = x.shape
    assert len(x) == len(y)
    for i in range(len(x)):
        assert x[i] == y[i]
    return True