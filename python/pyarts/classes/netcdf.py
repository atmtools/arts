#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 10:43:32 2020

@author: richard
"""


import pyarts
import netCDF4
import numpy as np


def check_data(record):
        name = record[0]
        data = record[1]
        typ = record[2]
        dims = record[3]
        
        if not isinstance(name, str):
            return False
        if not isinstance(typ, type):
            return False
        if typ == float:
            if not (isinstance(data, np.ndarray) or isinstance(data, float)):
                return False
        elif typ == str:
            if not isinstance(data, np.str):
                return False
        if not isinstance(dims, dict):
            return False
        return True    


def write_record_to_group(group, record):
    assert check_data(record), "Cannot understand output"
    name = record[0]
    data = record[1]
    typ = record[2]
    dims = record[3]
    
    curdims = []
    for dimname in dims:
        if dimname not in group.dimensions:
            group.createDimension(dimname, dims[dimname])
        curdims.append(dimname)
    
    if not len(curdims):
        if typ == str:
            group.setncattr_string(name, data)
        else:
            group.setncattr(name, data)
    else:
        var = group.createVariable(name, typ, curdims)
        var[:] = data


def write_to_group(group, ncdata):
    data, arraytype = ncdata
    
    if arraytype:
        group.nelem = len(data)
        for i in range(len(data)):
            innergroup = group.createGroup("pos{}".format(i))
            write_to_group(innergroup, data[i])
    else:
        for record in data:
            write_record_to_group(group, record)


def save_workspace_variable(nc, var):
    ncdata = pyarts.classes.from_workspace(var).netcdfify()

    group = nc.createGroup(var.name)
    group.type = var.group
    write_to_group(group, ncdata)
        

def save(filename, arts_list):
    nc = netCDF4.Dataset(filename, "w", format="NETCDF4")
    
    nc.pyarts_version = pyarts.version
    nc.creation_date = str(pyarts.classes.Time())
    try:
        for x in arts_list:
            save_workspace_variable(nc, x)
        nc.close()
    except:
        nc.close()
        raise RuntimeError("Cannot write all variables, closing...")


def load(filename):
    nc = netCDF4.Dataset("play.nc", "r")
    
    out = {}
    for key in nc.groups:
        group = nc[key]
        out[key] = pyarts.classes.__dict__[group.type]()
        out[key].denetcdf(group)
    return out
