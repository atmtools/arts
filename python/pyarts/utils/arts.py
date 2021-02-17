# -*- coding: utf-8 -*-

"""Collection of utility functions."""

import numpy as np

from pyarts.xml.names import basic_types, tensor_names, complex_tensor_names


def get_arts_typename(var):
    """Returns the ARTS type name for this variable.

    Args:
        var: Variable to get the ARTS type name for.

    Returns:
        str: ARTS type name.

    """
    if type(var).__name__ in basic_types:
        ret = basic_types[type(var).__name__]
        if ret == 'Array':
            if len(var) == 0:
                return None
            else:
                element_type = get_arts_typename(var[0])
                for element in var[1:]:
                    if element_type != get_arts_typename(element):
                        return None
                ret = 'ArrayOf' + element_type
    elif isinstance(var, np.ndarray):
        if np.issubdtype(var.dtype, np.complex128):
            ret = complex_tensor_names[var.ndim - 1]
        else:
            ret = tensor_names[var.ndim - 1]
    else:
        ret = type(var).__name__

    return ret


def return_if_arts_type(var, artstype):
    """If variable is of specified ARTS type, return its value.

    Parameters:
        var: variable to check
        artstype: wanted ARTS type

    Returns: value of var, if var is specified type. None for NoneType
    """
    arts_name = get_arts_typename(var)
    if arts_name is None:
        return None
    elif arts_name == artstype:
        return var
    else:
        raise TypeError('Expected {} but got {}.' .format(artstype, arts_name))
