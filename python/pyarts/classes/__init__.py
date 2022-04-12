# -*- coding: utf-8 -*-

"""This module contains classes within ARTS.
"""

# Full interface including internal manipulation
import pyarts.pyarts_cpp as cxx

for attr in dir(cxx):
    if isinstance(eval(f"cxx.{attr}"), type(cxx.Index)):
        exec(f"{attr} = cxx.{attr}")
del attr
