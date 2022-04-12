#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 10:43:32 2020

@author: richard
"""


import pyarts
import netCDF4
import numpy as np


def check_data(record: list) -> bool:
    """ Checks a netcdfifable record if it is good or not
    
    Throws if the record is not Sized:4
    
    Input:
        A list of a string, some data, the underlying data type, a dict of dims of ints
    
    Output:
        True or False
    """
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
    else:
        for key in dims:
            if not isinstance(dims[key], int):
                return False
    return True    


def write_record_to_group(group: netCDF4.Group, record: list) -> None:
    """ Writes a netcdfiable record to the group
    
    Input:
        group:
            A place to write to
        
        record:
            A list of a string, some data, the underlying data type, a dict of dims of ints
    """
    assert check_data(record), "Cannot understand netfifable record"
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


def write_to_group(group: netCDF4.Group, ncdata: list) -> None:
    """ Writes a list of netcdfiable records to the group or to array-like
    
    Input:
        group:
            A place to write to
        
        record:
            A tuple of lists and an array bool
    """
    data, arraytype = ncdata
    
    if arraytype:
        group.nelem = len(data)
        for i in range(len(data)):
            innergroup = group.createGroup("pos{}".format(i))
            write_to_group(innergroup, data[i])
    else:
        for record in data:
            write_record_to_group(group, record)


def save_workspace_variable(nc: netCDF4.Dataset, var: pyarts.workspace.WorkspaceVariable) -> None:
    """ Saves a variable to a netcdf dataset
    
    Input:
        nc:
            A writable netcdf dataset
        
        var:
            A workspace variable
    """
    ncdata = pyarts.classes.from_workspace(var).netcdfify()

    group = nc.createGroup(var.name)
    group.type = var.group
    write_to_group(group, ncdata)
        

def save(filename: str, arts_list: list) -> None:
    """ Saves a list of workspace variable to a netcdf file
    
    Tries to save the variable but close the file in a try-except catch
    
    Input:
        filename:
            A path, preferably a .nc file (note that .nc files blocks after bad writes)
        
        arts_list:
            A list of workspace variables
    
    Output:
        A file at path "filename"
    """
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


def load(filename: str) -> dict:
    """ Loads a netcdf file into a dictionary of pyarts.classes
    
    Input:
        filename:
            The name of a readable netcdf file
    
    Output:
        A dict of classes
    """
    nc = netCDF4.Dataset(filename, "r")
    
    out = {}
    for key in nc.groups:
        group = nc[key]
        out[key] = pyarts.classes.__dict__[group.type]()
        out[key].denetcdf(group)
    return out


def set_workspace(input: dict, arts: pyarts.workspace.Workspace) -> None:
    """ Set workspace from a dict
    
    All variables that cannot be read are added to the workspace as
    "newpython_KEY", where KEY is a dictionary key [uses staticmethod: name()]
    
    Input:
        input:
            A dictionary of pyarts.classes, where the dictionary name is a valid
            name of the same type in the ARTS workspace
    
        arts:
            A valid ARTS workspace
    """
    
    for key in input:
        try:
            exec("pyarts.classes.from_workspace(arts.{}).set(input[key])".format(key))  # FIXME:  there must be a better way to do this...
        except:
            print("Trying to create variable: newpython_{}".format(key))
            arts.create_variable(input[key].name(), "newpython_{}".format(key))
            exec("pyarts.classes.from_workspace(arts.newpython_{}).set(input[key])".format(key))  # ...
    