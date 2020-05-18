#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 10:43:32 2020

@author: richard
"""


import pyarts
import netCDF4
import numpy as np


def save_workspace_variable(nc, var):
    try:
        names, data, types, dims = pyarts.classes.from_workspace(var).netcdfify()
    except:
        raise ValueError("netcdify() must return correct types")
    
    group = nc.createGroup(var.name)
    
    group.type = var.group
    for i in range(len(names)):
        
        curdims = []
        for dimname in range(len(dims[i])):
            if dimname not in group.dimensions:
                group.createDimension(dimname, dims[i][dimname])
            curdims.append(dimname)
        
        group.createVariable(names[i], types[i], curdims)
        group.variables[names[i]][:] = data[i]


def save(filename, arts_list):
    nc = netCDF4.Dataset(filename, "w", format="NETCDF4")
    
    nc.pyarts_version = pyarts.version
    try:
        for x in arts_list:
            save_workspace_variable(nc, x)
        nc.close()
    except:
        nc.close()
        raise RuntimeError("Cannot write all variables, closing...")