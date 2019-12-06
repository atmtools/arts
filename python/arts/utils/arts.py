# -*- coding: utf-8 -*-

"""Collection of utility functions."""

import numpy as np

from arts.xml.names import basic_types, tensor_names, complex_tensor_names


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


def as_quantumnumbers(var):
    """Takes a quantum number prospect and turns it into a quantum number type
    if possible

    Parameters:
        var (dict, QuantumNumberRecord, QuantumNumbers, None, str):
            Quantum numbers

    Returns:
        QN (QuantumNumberRecord, QuantumNumbers): Returned quantumn numbers.
            No change if already quantum numbers type
    """

    if type(var) in [QuantumNumberRecord,
                     QuantumNumbers]:
        return var
    elif var is None:
        return QuantumNumbers()

    assert type(var) in [dict, str], "Cannot recognize as quantum number"

    if 'UP' in var and 'LO' in var:
        if type(var) is dict:
            return QuantumNumberRecord.from_dict(var)
        else:
            return QuantumNumberRecord.from_str(var)
    elif 'UP' in var:
        if type(var) is dict:
            var['LO'] = {}
            return QuantumNumberRecord.from_dict(var)
        else:
            return QuantumNumberRecord.from_str(var + ' LO')
    elif 'LO' in var:
        if type(var) is dict:
            var['UP'] = {}
            return QuantumNumberRecord.from_dict(var)
        else:
            return QuantumNumberRecord.from_str(var + ' UP')
    else:
        return QuantumNumbers(var)

from arts.catalogues import QuantumNumberRecord
from arts.catalogues import QuantumNumbers
