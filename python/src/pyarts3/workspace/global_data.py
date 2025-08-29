import numpy as np
import scipy as sp
import pyarts3.arts as cxx


def convert(group, value):
    """ Converts a value into something that can be turned into an Arts group
    
    This is intended to be used purely by a function parsing old controlfiles
    """
    groupname = str(group)
    
    if isinstance(value, str) or isinstance(value, cxx.String):
        try:
            eval(f"cxx.{groupname}('{value}')")
            # We understand the value already!
        except:
            value = eval(str(value))  # We try to stringify the value
    
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
