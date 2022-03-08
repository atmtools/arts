import numpy as np
import scipy as sp
import pyarts.pyarts_cpp as cxx

def get_raw_method_map():
    """ Returns a combined Arts MdRawMap and md_data_raw """
    map = cxx.get_MdRawMap()
    data = cxx.get_md_data_raw()

    for key in map:
        map[key] = data[map[key]]
    
    return map

def get_variables_map():
    """ Returns a combined Arts WsvMap and wsv_data """
    map = cxx.get_WsvMap()
    data = cxx.get_wsv_data()

    for key in map:
        map[key] = data[map[key]]
    
    return map

def get_variable_name(i):
    """
    Lookup the name of a variable given its workspace index.
    Args:
        i(int): The index of the workspace variable.
    Returns:
        str: The name of the workspace variable.
    """
    vars = cxx.get_wsv_data()
    return str(vars[i].name)

def convert(group, value):
    """ Converts a value into something that can be turned into an Arts group
    
    This is intended to be used purely by a function parsing old controlfiles
    """
    
    if isinstance(group, int):
        data = cxx.get_wsv_group_names()
        groupname = data[group]
    else:
        groupname = group
    groupname = str(groupname)
    
    if isinstance(value, cxx.String):
        value = str(value)
    
    if isinstance(value, str) and groupname != "String":
        value = eval(value)
    
    if groupname == "Index":
        return np.int64(value)
    
    if groupname == "ArrayOfIndex":
        return np.array(value, dtype=np.int64, order='C', ndmin=1)
    
    if groupname == "ArrayOfArrayOfIndex":
        return eval(np.array2string(np.array(value, dtype=np.int64, order='C', ndmin=2)))
        
    if groupname == "String":
        return str(value)
        
    if groupname == "Numeric":
        return np.float64(value)
    
    if groupname == "Vector":
        return np.array(value, dtype=np.float64, order='C', ndmin=1)
    
    if groupname == "Matrix":
        return np.array(value, dtype=np.float64, order='C', ndmin=2)
        
    if groupname == "Sparse":
        return sp.sparse.coo_matrix(value)
        
    if groupname.startswith("Tensor"):
        dim = int(group[6])
        return np.array(value, dtype=np.float64, order='C', ndmin=dim)
    
    return value
