#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 09:06:26 2023

@author: richard
"""


import pyarts
import numpy as np
import matplotlib.pyplot as plt


def all_keys(atm):
    return atm[0].keys() if len(atm) else []


def select_keys(atm, keys):
    return all_keys(atm) if len(keys) == 0 else keys


def select_atm(atm, keys=[]):
    keys = select_keys(atm, keys)

    return {key: np.array([x[key] for x in atm]) for key in keys} \
        if len(atm) else dict.fromkeys(keys)


class Info:
    def __init__(self, xscale, xlabel, xfactor) -> None:
        self.xscale = xscale
        self.xlabel = xlabel
        self.xfactor = xfactor


def plt_info(atm, keys=[]):
    keys = select_keys(atm, keys)
    
    out = {}
    for key in keys:
        if isinstance(key, pyarts.arts.ArrayOfSpeciesTag):
            out[key] = Info("linear", f"VMR {key} [-]", 1)
        elif isinstance(key, pyarts.arts.options.AtmKey):
            if key == pyarts.arts.options.AtmKey.t:
                out[key] = Info("linear", "Temperature [K]", 1)
            elif key == pyarts.arts.options.AtmKey.p:
                out[key] = Info("log", "Pressure [Pa]", 1)
            elif key == pyarts.arts.options.AtmKey.mag_u:
                out[key] = Info("linear", "Magnetic u-component [µT]", 1e6)
            elif key == pyarts.arts.options.AtmKey.mag_v:
                out[key] = Info("linear", "Magnetic v-component [µT]", 1e6)
            elif key == pyarts.arts.options.AtmKey.mag_w:
                out[key] = Info("linear", "Magnetic w-component [µT]", 1e6)
            elif key == pyarts.arts.options.AtmKey.wind_u:
                out[key] = Info("linear", "Wind u-component [m/s]", 1)
            elif key == pyarts.arts.options.AtmKey.wind_v:
                out[key] = Info("linear", "Wind v-component [m/s]", 1)
            elif key == pyarts.arts.options.AtmKey.wind_w:
                out[key] = Info("linear", "Wind w-component [m/s]", 1)
            else:
                assert False, "Unknown key type"
        else:
            assert False, "Unknown key type"
    return out


def ppvar_atm(atm, path, keys=[], info={}):
    y = 1.0 * path.pos[:, 0]

    xs = select_atm(atm, keys)

    info = info if len(info) else plt_info(atm, keys)

    for key in xs:
        s = info[key]
        x = xs[key]
        
        plt.plot(s.xfactor * x, y)
        plt.xscale(s.xscale)
        plt.xlabel(s.xlabel)
        plt.ylabel("Altitude [m]")
        
        plt.show()
