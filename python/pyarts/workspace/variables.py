""" The variables submodule.

This module contains symbolic representations of all ARTS workspace variables.

The variables are loaded dynamically when the module is imported, which ensures that they
up to date with the current ARTS build.

TODO: The group names list is redudant w.rt. group_ids.keys(). Should be removed.

Attributes:
    group_names([str]): List of strings holding the groups of ARTS WSV variables.
    group_ids(dict):    Dictionary mapping group names to the group IDs which identify
                        groups in the ARTS C API.
"""

import ctypes as c
import os
import numpy as np
import re
import scipy as sp
import tempfile
import weakref

from pyarts.workspace.api import arts_api
from pyarts.workspace.agendas import Agenda
from pyarts.xml.names import tensor_names


class WorkspaceVariable:
    """
    The WorkspaceVariable represents ARTS workspace variables in a symbolic way. This
    means that they are not associated with a single workspace and therefore do not have a
    unique value. Their value in a given workspacecan be accessed, however, using the value()
    method.

    Attributes:
        ws_id(int): The Index variable identifying the variable in the ARTS C API.
        name(str): The name of the workspace variable.
        group(str): The name of the group this variable belongs to.
        description(str): The documentation of the variable as in methods.cc
        ws(Workspace): The associated workspace (if available) or None
    """
    def __init__(self, ws_id, name, group, description, ws = None):
        self.ws_id = ws_id
        self.name = name
        self.group = group
        self.group_id = group_ids[group]
        self.description = description
        if not ws is None:
            self._ws = weakref.ref(ws)
        else:
            self._ws = None

        self.ndim = None
        if self.group == "Vector":
            self.ndim = 1
        if self.group == "Matrix":
            self.ndim = 2
        m = re.match(r"^Tensor(\d)$", self.group)
        if m:
            self.ndim = int(m.group(1))

        self.update()

    @property
    def ws(self):
        if self._ws is None:
            return None
        else:
            return self._ws()

    @ws.setter
    def ws(self, workspace):
        if workspace is None:
            self._ws = workspace
        else:
            self._ws = weakref.ref(workspace)

    def __getstate__(self):
        return self.ws_id, self.name, self.group, \
               self.group_id, self.description, self.ndim

    def __setstate__(self, state):
        self.ws_id, self.name, self.group, self.group_id, self.description,\
            self.ndim = state

    def __repr__(self):
        return "ARTS Workspace Variable: " + self.name

    def help(self):
        """
        Print ARTS documentation for this variable.
        """
        s  = "ARTS Workspace Variable\n\n"
        s += "Name:  " + self.name + "\n"
        s += "Group: " + self.group + "\n\n"
        s += self.description
        print(s)
        return None

    def __str__(self):
        return self.__repr__()

    def __setattr__(self, name, value):

        if name == "value":
            if self.ws is None:
                raise Exception("Cannot set value of WSV without associated "
                                " workspace.")
            else:
                self.ws.__setattr__(self.name, value)
        else:
            super().__setattr__(name, value)

    def print(self):
        """ Print variable value using ARTS Print(...) WSM.

        Raises:
            Exception: If the variable has no associated workspace.
        """
        if (self.ws):
            self.ws.Print(self, 0)
        else:
            raise Exception("Can't print variable without associated ARTS workspace.")

    @staticmethod
    def get_variable_name(i):
        """
        Lookup the name of a variable given its workspace index.

        Args:
            i(int): The index of the workspace variable.

        Returns:
            str: The name of the workspace variable.
        """
        s = arts_api.get_variable(i)
        name = s.name.decode("utf8")
        return name

    @staticmethod
    def get_group_id(value):
        """ This static method is used to determine how (and if) a given python variable can
        be mapped to a ARTS workspace variable group. The returned group id is required to
        add the variable to a workspace.

        Args:
            value(any): The python variable to map to the ARTS group.

        Returns:
            int: The index of the group which can be used to represent the python variable
                 or None if the type is not supported.
        """
        import numpy as np
        if isinstance(value, WorkspaceVariable):
            return group_ids[value.group]
        elif isinstance(value, Agenda):
            return group_ids["Agenda"]
        elif isinstance(value, np.long):
            return group_ids["Index"]
        elif isinstance(value, (float, np.double)):
            return group_ids["Numeric"]
        elif isinstance(value, str):
            return group_ids["String"]
        elif isinstance(value, np.ndarray):
            if value.ndim == 1:
                return group_ids["Vector"]
            elif value.ndim == 2:
                return group_ids["Matrix"]
            elif value.ndim == 3:
                return group_ids["Tensor3"]
            elif value.ndim == 4:
                return group_ids["Tensor4"]
            elif value.ndim == 5:
                return group_ids["Tensor5"]
            elif value.ndim == 6:
                return group_ids["Tensor6"]
            elif value.ndim == 7:
                return group_ids["Tensor7"]
            else:
                raise ValueError(
                    "Numpy arrays are only supported up to 7 dimensions."
                )
        elif sp.sparse.issparse(value):
            return group_ids["Sparse"]
        elif type(value) == list:
            group_name = ""
            nested_value = value
            while type(nested_value) == list and len(nested_value) > 0:
                nested_value = nested_value[0]
                group_name += "ArrayOf"
            if type(nested_value) == list and len(nested_value) == 0:
                raise ValueError("Empty lists are currently not handled.")
            else:
                typename = type(nested_value).__name__
                if isinstance(nested_value, str):
                    group_name += "String"
                    return group_ids[group_name]
                elif isinstance(nested_value, np.long):
                    group_name += "Index"
                    return group_ids[group_name]
                elif isinstance(nested_value, (float, np.double)):
                    raise ValueError("Vectors, Matrices or Tensors should be"
                                     " passed as numpy.ndarray and not as"
                                     " lists.")
                elif hasattr(nested_value, 'write_xml') and typename in group_names:
                    return group_ids[group_name + typename]
                elif isinstance(nested_value, np.ndarray):
                    group_name += tensor_names[len(nested_value.shape) - 1]
                    return group_ids[group_name]
                else:
                    raise ValueError(
                        f"Nested array with internal type "
                        f"{type(nested_value)} not supported.")
        elif hasattr(value, 'write_xml') and type(value).__name__ in group_names:
            return group_ids[type(value).__name__]
        else:
            raise ValueError(f"Type {type(value)} currently not supported.")

    @classmethod
    def convert(cls, group, value):
        """ Tries to convert a given python object to an object of the python class
        representing the given ARTS WSV group.

        Args:
            group(string): The name of an ARTS WSV group.
            group(any):    The object to convert

        Returns:
            (any): The converted object or None is conversion was unsuccessful.
        """
        import numpy as np

        try:
            gid = cls.get_group_id(value)
            if group_names[gid] == group:
                return value
        except:
            pass

        if (group == "Index"):
            return int(value)
        if (group == "String"):
            return str(value)
        if (group == "ArrayOfString"):
            return [str(i) for i in value]
        if (group == "Numeric"):
            return np.float64(value)
        if (group == "Vector"):
            return np.array(value, dtype=np.float64, order='C', ndmin=1)
        if (group == "Matrix"):
            return np.array(value, dtype=np.float64, order='C', ndmin=2)
        if (group == "Sparse"):
            return sp.sparse.coo_matrix(value)
        if (group[:6] == "Tensor"):
            dim = int(group[6])
            return np.array(value, dtype=np.float64, order='C', ndmin=dim)
        if group.startswith("ArrayOf"):
            subgroup = group[7:]
            if hasattr(value, "__iter__"):
                return [cls.convert(subgroup, v) for v in value]
            else:
                return [cls.convert(subgroup, value)]
        return None

    @staticmethod
    def iter():
        """
        Iterator returning a WorkspaceVariable object for each ARTS WSV available.
        """
        for i in range(arts_api.get_number_of_variables()):
            s = arts_api.get_variable(i)
            name        = s.name.decode("utf8")
            description = s.description.decode("utf")
            group       = group_names[s.group]
            yield WorkspaceVariable(i, name, group, description)

    @property
    def initialized(self):

        ws = self.ws
        if ws is None:
            raise ValueError("WorkspaceVariable object needs associated"
                             " Workspace to determine value.")

        v = arts_api.get_variable_value(ws.ptr, self.ws_id, self.group_id)
        return v.initialized

    @property
    def value(self):
        """ Return the value of the variable in a given workspace.

        By default this function will check the value in the workspace associated
        with the variable of in the workspace object provided as argument to the
        function call. If the variable has an associated workspace the workspace
        provided as argument will be ignored.

        Returns:
            The value of the workspace variable represented by an object of
            the corresponding python types.

        Raises:
            Exception: If the type of the workspace variable is not supported
            by the interface.

        """
        from pyarts.types import classes as arts_classes


        if (self.ws):
            ws = self.ws
        if not ws:
            raise ValueError("WorkspaceVariable object need Workspace to determine value.")

        v = arts_api.get_variable_value(ws.ptr, self.ws_id, self.group_id)
        if not v.initialized:
            raise Exception("WorkspaceVariable " + self.name + " is uninitialized.")

        if self.group in arts_classes:
            cls = arts_classes[self.group]
            if hasattr(cls, "__from_variable_value_struct__"):
                return cls.__from_variable_value_struct__(v)
        if self.group == "Index":
            return c.cast(v.ptr, c.POINTER(c.c_long))[0]
        elif self.group == "Numeric":
            return c.cast(v.ptr, c.POINTER(c.c_double))[0]
        elif self.group == "String":
            return (c.cast(v.ptr, c.c_char_p)).value.decode("utf8")
        elif self.group == "ArrayOfIndex":
            return [c.cast(v.ptr, c.POINTER(c.c_long))[i]
                    for i in range(v.dimensions[0])]
        elif self.group == "Sparse":
            m    = v.dimensions[0]
            n    = v.dimensions[1]
            nnz  = v.dimensions[2]
            if nnz == 0:
                return sp.sparse.csr_matrix(0)
            else:
                data = np.ctypeslib.as_array(c.cast(v.ptr,
                                                    c.POINTER(c.c_double)),
                                             (nnz,))
                row_indices = np.ctypeslib.as_array(v.inner_ptr, (nnz,))
                col_starts  = np.ctypeslib.as_array(v.outer_ptr, (m + 1,))
                return sp.sparse.csr_matrix((data, row_indices, col_starts),
                                            shape=(m,n))
        elif self.group == "Agenda":
            return Agenda(v.ptr)
        elif self.ndim:
            shape = []
            size  = 1
            for i in range(self.ndim):
                shape.append(v.dimensions[i])
                size *= v.dimensions[i]
            if size > 0:
                self.__array_interface__ = {"shape"  : tuple(shape),
                                            "typestr" : "|f8",
                                            "data" : (v.ptr, False),
                                            "version" : 3}
                return np.asarray(self)
            else:
                return np.zeros(shape)
        else:
            try:
                return self.to_arts()
            except:
                raise Exception("Type of workspace variable is not supported "
                                + " by the interface.")

    def update(self):
        """ Update data references of the object.

        References to vector, matrices and tensors may change and must therefore
        be updated dynamically to ensure they are consistent with the state of
        the associated workspace. This method takes care of that.

        """
        if not self.ws==None and self.ndim:
            v = arts_api.get_variable_value(self.ws.ptr, self.ws_id, self.group_id)
            shape = []
            for i in range(self.ndim):
                shape.append(v.dimensions[i])
            self.__array_interface__ = {"shape"  : tuple(shape),
                                        "typestr" : "|f8",
                                        "data" : (v.ptr, False),
                                        "version" : 3}

    def erase(self):
        """
        Erase workspace variable from its associated workspace.
        """
        if self.ws:
            arts_api.erase_variable(self.ws.ptr, self.ws_id, self.group_id)
            self.ws = None

    def describe(self):
        """
        Print the description of the variable as given in ARTS methods.cc
        """
        print(self.description.format())

    def to_arts(self):
        """
        Return the value of this variable as a typhon type. This function
        writes the value of the variable to a temporary file and reads it
        into Python using typhon load function. The purpose of this function
        is to access WSV whose groups are not natively supported by the
        C API.

        Returns:
            A typhon object with the same value as the WSV in the associated
            workspace.
        """
        from pyarts.xml import load
        from pyarts.types import classes as arts_classes
        from pyarts import classes as native_classes

        if not self.ws:
            raise Exception("Cannot retrieve the value of a variable without "
                            + " associated Workspace.")
        with tempfile.TemporaryDirectory() as tmpdir:
            tfile = os.path.join(tmpdir, 'wsv.xml')
            self.ws.WriteXML("binary", self, tfile)
            if self.group in arts_classes:
                v = load(tfile)
            else:
                cls = getattr(native_classes, self.group)
                v = cls()
                v.readxml(tfile)

        return v

    def from_arts(self, var):
        """
        Set the value of this WSV in the associated workspace to the given
        typhon type. This function writes the value in ASCII format to a
        temporary file and reads it into the workspace

        Args:
            var: The value to which this WSV should be set in the associated
                 workspace.

        """
        from pyarts.xml import save

        if not self.ws:
            raise Exception("Cannot set the value of a variable without "
                            + " associated Workspace.")
        with tempfile.TemporaryDirectory() as tmpdir:
            tfile = os.path.join(tmpdir, 'wsv.xml')
            save(var, tfile, format='binary')
            self.ws.ReadXML(self, tfile)


# Get ARTS WSV groups
group_names = [arts_api.get_group_name(i).decode("utf8")
               for i in range(arts_api.get_number_of_groups())]
group_ids = dict([(id, name) for (name,id) in enumerate(group_names)])


workspace_variables = dict()
for v in WorkspaceVariable.iter():
    workspace_variables[v.name] = v

