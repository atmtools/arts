#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 09:06:26 2023

@author: richard
"""


import numpy as np
import matplotlib.pyplot as plt


def select_atm(atm, keys=[]):
    keys = atm[0].keys() if len(atm) and len(keys) == 0 else keys

    return {key: np.array([x[key] for x in atm]) for key in keys} \
        if len(atm) else dict.fromkeys(keys)


def ppvar_atm(atm, path, var=[]):
    pos = 1.0 * path.pos[:, 0]
    x = select_atm(atm, var)

    if not isinstance(x, dict):
        x = {var: x}

    for key in x:
        plt.plot(x[key], pos)
        plt.show()
