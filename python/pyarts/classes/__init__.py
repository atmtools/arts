# -*- coding: utf-8 -*-

"""This module contains classes within ARTS.
"""

# Attempt at conversions
from pyarts.workspace.api import arts_api as lib, c

from pyarts_cpp.classes import *  #noqa

def from_workspace(x):
    """ Converts a workspace variable to a python class

    Input:
        x:
            A workspace variable (pyarts.workspace.WorkspaceVariable)

    Output:
        A modifyable ARTS variable (pyarts.workspace.eval(x.group))
    """
    import pyarts
    import pyarts_cpp
    if isinstance(x, pyarts.workspace.WorkspaceVariable):
        v = lib.get_variable_data_pointer(x.ws.ptr, x.ws_id)
        
        try:
            retfun = eval(f"pyarts_cpp.workspace_references.{x.group}")
        except:
            raise RuntimeError(f"Cannot find reference for type {x.group}")
        
        return retfun(v)
    else:
        raise TypeError("Expects WorkspaceVariable")


lib.get_variable_data_pointer.restype = c.c_void_p
lib.get_variable_data_pointer.argtypes = [c.c_void_p, c.c_long]
