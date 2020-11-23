"""ARTS C API Interface

This module provides a foreign function interface for the ARTS C API.
It defines the C structs used by the interface as ctypes.Structure
child classes as well as the return argument and return types of the
function provided by the C API.

Requirements
------------

The ARTS C API is provided by the arts_api.so library and is required by
the module. The module will check if the ``ARTS_BUILD_PATH`` variable is set
and assume the library can be found in the src subdirectory. If opening
the library fails loading the module will fail with an EnvironmentError.

Attributes:

    arts_api(CDLL): The ctypes library handle holding the ARTS C API.

"""

import ctypes as c
import logging
import os
import pkg_resources

import numpy as np
import scipy as sp
import sys

from pyarts.environment import environ


logger = logging.getLogger(__name__)

################################################################################
# Version Requirements
################################################################################

arts_minimum_major    = 2
arts_minimum_minor    = 3
arts_minimum_revision = 1167

################################################################################
# Load ARTS C API
################################################################################

try:
    lib_path = pkg_resources.resource_filename(__name__, "libarts_api.so")
    arts_api = c.cdll.LoadLibrary(lib_path)
except Exception as e:
    raise EnvironmentError("The following error was encountered when "
                           "trying to load the ARTS API: ", e)

################################################################################
# Version Check
################################################################################

class VersionStruct(c.Structure):
    """
    The ARTS version is represented by 3 values of type long: the major, minor
    and revision number.
    """
    _fields_ = [("major", c.c_long),
                ("minor", c.c_long),
                ("revision", c.c_long)]

arts_api.get_version.argtypes = None
arts_api.get_version.restype  = VersionStruct

version = arts_api.get_version()
if (version.major, version.minor, version.revision) \
   < (arts_minimum_major, arts_minimum_minor, arts_minimum_revision):

    raise EnvironmentError("This arts version requires at least arts-"
                           + str(arts_minimum_major) + "."
                           + str(arts_minimum_minor) + "."
                           + str(arts_minimum_revision) + " of ARTS.")

################################################################################
# Initialize API
################################################################################

arts_api.initialize()
try:
    filename = sys.modules["__main__"].__file__
    basename, _ = os.path.splitext(os.path.basename(filename))
    arts_api.set_basename(c.c_char_p(basename.encode()))
except:
    pass

################################################################################
# ARTS runtime environment manipulation
################################################################################

def find_controlfile(name):
    """ Recursively search arts include path for given file.
    Args:
        name(str): Name of the file.
    Raises:
        Exception: If the file cannot be found.
    Returns:
        path(str): The full path of the file.
    """
    paths = arts_include_path + [os.getcwd()]
    path  = None

    for p in paths:
        if os.path.isfile(os.path.join(p, name)):
            path = os.path.join(p, name)
    if (path):
        return path
    else:
        raise Exception("File " + name + " not found. Search path was:\n "
                        + str(paths))

def include_path_push(path):
    """
    Add path to include path of the ARTS runtime.

    Args:
        path(str): Path to add to the ARTS include path.
    """
    arts_api.include_path_push(c.c_char_p(path.encode()))

def include_path_pop():
    """
    Remove most recently added include path.
    """
    arts_api.include_path_pop()

def data_path_push(path):
    """
    Add path to data path of the ARTS runtime.

    Args:
        path(str): Path to add to the ARTS data path.
    """
    arts_api.data_path_push(c.c_char_p(path.encode()))

def data_path_pop():
    """
    Remove most recently added data path.
    """
    arts_api.data_path_pop()

################################################################################
# Python values that represent empty elements.
################################################################################

def is_empty(value):
    return (value is None) or ((type(value) is list) and (len(value) == 0))

################################################################################
# ctypes Structures
################################################################################

class VariableStruct(c.Structure):
    """
    A (symbolic) ARTS workspace variable is represented using a struct containing
    pointers to the name and description of the method as well as the group id,
    i.e. the Index variable encoding the type of the variable.
    """
    _fields_ = [("name", c.c_char_p),
                ("description", c.c_char_p),
                ("group", c.c_long)]


class VariableValueStruct(c.Structure):
    """
    The ARTS C API uses C-structs to transfer workspace variables from and to the
    ARTS runtime. The VariableValueStruct class is used to access and create these
    structs in Python code. The fields of the C-struct can be accessed directly as
    attributes of the VariableValueStruct object.
    """
    _fields_ = [("ptr", c.c_void_p),
                ("initialized", c.c_bool),
                ("dimensions", 7 * c.c_long),
                (("inner_ptr"), c.POINTER(c.c_int)),
                (("outer_ptr"), c.POINTER(c.c_int))]

    @classmethod
    def empty(cls):
        s = VariableValueStruct(0)
        s.ptr = 0
        return s

    def __init__(self, value):
        """ Create a VariableValue struct from a python object.

        This functions creates a variable value struct from a python object so
        that it can be passed to the C API. If the type of the object is not
        supported, the data pointer will be NONE.

        The built-in Python types that are currently supported are:

            - int
            - float
            - string
            - lists of int and lists of string
            - numpy.ndarray
            - scipy.sparse

        User defined classes are supported through a generic interface. The constructor
        looks for an attribute function __to_value_struct__, which should return a dictionary
        containing the value associated with the fields of the C-struct.

        Args:
            value(object): The python object to represent as a VariableValue struct.

        """
        ptr = None
        initialized = True
        dimensions  = [0] * 7

        self._temporaries = []

        # Generic interface
        if hasattr(value, "__to_value_struct__"):
            d = value.__to_value_struct__()
            if "ptr" in d:
                ptr = d["ptr"]
            if "dimensions" in d:
                dimensions = d["dimensions"]
        # Index
        elif isinstance(value, np.long):
            ptr = c.cast(c.pointer(c.c_long(value)), c.c_void_p)
        # Numeric
        elif isinstance(value, (float, np.double)):
            temp = np.float64(value)
            ptr = c.cast(c.pointer(c.c_double(temp)), c.c_void_p)
        # String
        elif isinstance(value, str):
            ptr = c.cast(c.c_char_p(value.encode()), c.c_void_p)
        # Vector, Matrix
        elif isinstance(value, np.ndarray):
            # arrays need to be contiguous when passed to the ARTS API
            value = np.ascontiguousarray(value)

            if value.dtype == np.float64:
                ptr = value.ctypes.data
                for i in range(value.ndim):
                    dimensions[i] = value.shape[i]
        # Scipy sparse matrices
        elif sp.sparse.issparse(value):
            m = sp.sparse.coo_matrix(value)
            self._temporaries += [m]
            dimensions[0] = m.shape[0]
            dimensions[1] = m.shape[1]
            dimensions[2] = m.nnz
            ptr = m.data.ctypes.data
            self.inner_ptr = c.cast(m.row.ctypes.data, c.POINTER(c.c_int))
            self.outer_ptr = c.cast(m.col.ctypes.data, c.POINTER(c.c_int))
        # Array of String or Integer
        elif isinstance(value, list):
            if not value:
                raise ValueError("Empty lists currently not supported.")
            ps = []
            if isinstance(value[0], str):
                for s in value:
                    ps.append(c.cast(c.c_char_p(s.encode()), c.c_void_p))
                p_array = (c.c_void_p * len(value))(*ps)
                ptr = c.cast(c.pointer(p_array), c.c_void_p)
            if isinstance(value[0], np.long):
                ptr = c.cast(c.pointer((c.c_long * len(value))(*value)), c.c_void_p)
            dimensions[0] = len(value)

        self._value = value
        self.ptr = ptr
        self.initialized = initialized
        self.dimensions  = (c.c_long * 7)(*dimensions)

class CovarianceMatrixBlockStruct(c.Structure):
    """
    c struct representing block of covariance matrices.

    The indices field holds the indices of the corresponding
    retrieval quantities.

    The position field holds the row and column indices of the
    left- and upper-most element of the block w.r.t. the full
    covariance matrix.

    The dimension field hold the number of rows and columns of
    the block.

    The ptr field hold the pointer to the dense matrix data or
    to the element pointer of the sparse matrix that represents
    the block.

    The inner and outer pointer fields are null if the block is
    represented by a dense matrix. Otherwise these contain the
    pointers to the index array of the sparse matrix of which the
    block consists.
    """
    _fields_ = [("indices", 2 * c.c_long),
                ("position", 2 * c.c_long),
                ("dimensions", 2 * c.c_long),
                ("ptr", c.c_void_p),
                ("nnz", c.c_long),
                ("inner_ptr", c.POINTER(c.c_int)),
                ("outer_ptr", c.POINTER(c.c_int))]

class MethodStruct(c.Structure):
    """
    The method struct holds the internal index of the method (id), pointers
    to the null-terminated strings holding name and description, the number
    of generic inputs (n_g_in) and a pointer to the array holding the group ids
    of the output types, as well as the number of generic outputs and their types.
    """
    _fields_ = [("id", c.c_ulong),
                ("name", c.c_char_p),
                ("description", c.c_char_p),
                # Output
                ("n_out", c.c_ulong),
                ("outs", c.POINTER(c.c_long)),
                # Generic Output
                ("n_g_out", c.c_ulong),
                ("g_out_types", c.POINTER(c.c_long)),
                # Input
                ("n_in", c.c_ulong),
                ("ins", c.POINTER(c.c_long)),
                # Generic Input
                ("n_g_in", c.c_ulong),
                ("g_in_types", c.POINTER(c.c_long))]

################################################################################
# Function Arguments and Return Types
################################################################################

# Create ArtsWorkspace and return handle.
arts_api.create_workspace.argtypes = [c.c_long, c.c_long]
arts_api.create_workspace.restype  = c.c_void_p

# Destroy ArtsWorkspace instance from handle.
arts_api.destroy_workspace.argtypes = [c.c_void_p]
arts_api.destroy_workspace.restype   = None

# Include path manipulation.
arts_api.include_path_push.restype = None
arts_api.include_path_push.argtypes = [c.c_char_p]

# Data path manipulation.
arts_api.include_path_pop.restype = None
arts_api.include_path_pop.argtypes = None

arts_api.data_path_push.restype = None
arts_api.data_path_push.argtypes = [c.c_char_p]

arts_api.data_path_pop.restype = None
arts_api.data_path_pop.argtypes = None

# Set include ad data path of the arts runtime.
arts_api.get_error.restype  = c.c_char_p
arts_api.get_error.argtypes = None

arts_api.set_basename.restype  = None
arts_api.set_basename.argtypes = [c.c_char_p]

# Agendas
#
#
arts_api.create_agenda.argtypes = [c.c_char_p]
arts_api.create_agenda.restype  = c.c_void_p

arts_api.agenda_add_method.argtypes = [c.c_void_p, c.c_long,
                                       c.c_ulong, c.POINTER(c.c_long),
                                       c.c_ulong,c.POINTER(c.c_long)]
arts_api.agenda_add_method.restype  = None

arts_api.agenda_clear.argtypes = [c.c_void_p]
arts_api.agenda_clear.restype  = None

arts_api.agenda_insert_set.argtypes = [c.c_void_p, c.c_void_p, c.c_long, c.c_long]
arts_api.agenda_insert_set.restype = None

arts_api.agenda_append.argtypes = [c.c_void_p, c.c_void_p]
arts_api.agenda_append.restype = None

arts_api.parse_agenda.argtypes = [c.c_char_p]
arts_api.parse_agenda.restype  = c.c_void_p

arts_api.execute_agenda.argtypes = [c.c_void_p, c.c_void_p]
arts_api.execute_agenda.restype  = c.c_char_p

arts_api.destroy_agenda.argtypes = [c.c_void_p]
arts_api.destroy_agenda.restype  = None

# Groups
#
# Returns the number of WSV groups.
arts_api.get_number_of_groups.argtypes = None
arts_api.get_number_of_groups.restype  = c.c_ulong

# Return pointer to the name of the group with given index.
arts_api.get_group_name.argtypes = [c.c_long]
arts_api.get_group_name.restype  = c.c_char_p

# Variables
#
# Returns the number of (symbolic) workspace variable.
arts_api.get_number_of_variables.restype  = c.c_ulong
arts_api.get_number_of_variables.argtypes = None

# Returns workspace variable with index c_long as VariableStruct.
arts_api.lookup_workspace_variable.argtypes = [c.c_char_p]
arts_api.lookup_workspace_variable.restype  = c.c_long

# Returns workspace variable with index c_long as VariableStruct.
arts_api.get_variable.argtypes = [c.c_long]
arts_api.get_variable.restype  = VariableStruct

# Return pointer to variable value in a given workspace in the form of a VariableValueStruct.
arts_api.get_variable_value.argtypes = [c.c_void_p, c.c_long, c.c_long]
arts_api.get_variable_value.restype  = VariableValueStruct

# Set variable value in workspace given a workspace handle, the variable id, the group id
# and a VariableValueStruct
arts_api.set_variable_value.argtypes = [c.c_void_p, c.c_long, c.c_long, VariableValueStruct]
arts_api.set_variable_value.restype  =  c.c_char_p

# Adds a value of a given group to a given workspace.
arts_api.add_variable.restype  = c.c_long
arts_api.add_variable.argtypes = [c.c_void_p, c.c_long, c.c_char_p]

# Remove given variable from workspace.
arts_api.erase_variable.restype  = None
arts_api.erase_variable.argtypes = [c.c_void_p, c.c_long, c.c_long]

# Methods
#
# Returns the number of (symbolic) workspace variable.
arts_api.get_number_of_methods.restype  = c.c_ulong
arts_api.get_number_of_methods.argtypes = None

# Returns workspace variable with index c_long as VariableStruct.
arts_api.get_method.argtypes = [c.c_long]
arts_api.get_method.restype  = MethodStruct

# Return Pointer to name of jth generic output parameter of a given WSM.
arts_api.get_method_g_out.argtypes = [c.c_long, c.c_long]
arts_api.get_method_g_out.restype  = c.c_char_p

# Return Pointer to name of jth generic input parameter of a given WSM.
arts_api.get_method_g_in.argtypes = [c.c_long, c.c_long]
arts_api.get_method_g_in.restype  = c.c_char_p

# Return pointer to the default value of the jth generic input of a given WSM.
arts_api.get_method_g_in_default.argtypes = [c.c_long, c.c_long]
arts_api.get_method_g_in_default.restype  = c.c_char_p

# Return pointer to the default value of the jth generic input of a given WSM.
arts_api.get_g_in_nodef.argtypes = []
arts_api.get_g_in_nodef.restype  = c.c_char_p

# Return block from covariance matrix.
arts_api.get_covariance_matrix_block.argtypes = [c.c_void_p, c.c_long, c.c_bool]
arts_api.get_covariance_matrix_block.restype = CovarianceMatrixBlockStruct

# Execute a given workspace method.
arts_api.execute_workspace_method.argtypes = [c.c_void_p,
                                              c.c_long,
                                              c.c_ulong,
                                              c.POINTER(c.c_long),
                                              c.c_ulong,
                                              c.POINTER(c.c_long)]
arts_api.execute_workspace_method.restype  = c.c_char_p

# Print method documentation.
arts_api.method_print_doc.argtypes = [c.c_long]
arts_api.method_print_doc.restype  = c.c_char_p

# Callback insertion
arts_api.callbacks = []
arts_api.agenda_insert_callback.argtypes = [c.c_void_p, c.c_void_p]
arts_api.agenda_insert_callback.restype  = None

################################################################################
# Setup ARTS Environment
################################################################################

try:
    arts_include_path = environ.get("ARTS_INCLUDE_PATH").split(":")
except:
    arts_include_path = []
arts_include_path += [os.path.join("@ARTS_SRC_DIR@", "controlfiles")]

basedir = os.path.join(os.path.dirname(__file__), "..")
arts_include_path += [os.path.join(basedir, "controlfiles")]

try:
    arts_data_path = environ.get("ARTS_DATA_PATH").split(":")
except:
    arts_data_path = []
arts_data_path += [@ARTS_XML_DIR@]

# Set runtime parameters
for p in arts_include_path:
    include_path_push(p)

for p in arts_data_path:
    data_path_push(p)

include_path_push(os.getcwd())
data_path_push(os.getcwd())

nodef = arts_api.get_g_in_nodef().decode("utf8")
