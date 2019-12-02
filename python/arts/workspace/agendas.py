"""The agendas submodule.

This module provides the Agenda class, which is used to represent parsed
controlfiles and can be executed on a given `Workspace` object.

"""

import ctypes as c
import numpy  as np
import os

from arts.workspace.api import find_controlfile, arts_api
from arts.workspace.output import CoutCapture

class Agenda:
    def __init__(self, ptr):
        """ Initialize Agenda object from pointer to C API Agenda object.
        Args:
            ptr(c.c_void_p): Pointer to Agenda object created with the parse_agenda
            method of the ARTS C API.
        """
        self.ptr       = ptr

    @classmethod
    def create(cls, name):
        """
        Create agenda with given name.

        Parameters:

            name(str): The name of the agenda to create.

        Returns:

            The newly created agenda object with the given name.

        """
        ptr = arts_api.create_agenda(name.encode())
        return Agenda(ptr)

    def clear(self):
        """ Reset agenda to be empty."""
        arts_api.agenda_clear(self.ptr)

    def append(self, agenda):
        """Append agenda to this agenda.

        Parameters:

            agenda(:class:`Agenda`): The agenda to append to this agenda.

        Raises:

            Exception if :code:`agenda` is not of the :class:`Agenda`
        """

        if not isinstance(agenda, Agenda):
            raise Exception("Agenda to append must be of type Agenda.")

        arts_api.agenda_append(self.ptr, agenda.ptr)

    def add_method(*args, **kwargs):
        """
        Add a workspace method call to the agenda.

        Parameters:

            ws(arts.workspace.Workspace): A (dummy) workspace object.

            m(arts.workspace.WorkspaceMethod): The method to add to the
            agenda

            *args:  Positional arguments of the WSM call to add to the agenda.

            **kwargs: Key-word arguments of the WSM call to add to the agenda.
        """
        from arts.workspace.variables import group_names

        if len(args) < 3:
            raise Exception("Need at least self, a workspace and the method to"
                            " add as arguments.")
        self = args[0]
        ws   = args[1]
        m    = args[2]
        m_id, args_out, args_in, temps = m._parse_output_input_lists(ws,
                                                                     args[3:],
                                                                     kwargs)
        arg_out_ptr = c.cast((c.c_long * len(args_out))(*args_out),
                             c.POINTER(c.c_long))
        arg_in_ptr = c.cast((c.c_long * len(args_in))(*args_in),
                            c.POINTER(c.c_long))
        if not m.name[-3:] == "Set" or not m.name[:-3] in group_names:

            for t in temps:
                arts_api.agenda_insert_set(ws.ptr, self.ptr, t.ws_id, t.group_id)

            arts_api.agenda_add_method(c.c_void_p(self.ptr), m_id,
                                       len(args_out), arg_out_ptr,
                                       len(args_in), arg_in_ptr)
        else:
            from arts.workspace.variables import WorkspaceVariable

            name_out = WorkspaceVariable.get_variable_name(args_out[0])
            name_in = WorkspaceVariable.get_variable_name(args_in[0])
            wsv_out = getattr(ws, name_out)
            wsv_in  = getattr(ws, name_in)

            ws.Copy(wsv_out, wsv_in)

            group_id = arts_api.get_variable(args_out[0]).group
            arts_api.agenda_insert_set(ws.ptr, self.ptr, args_out[0], group_id)

    def add_callback(self, f):
        """
        Add a Python callback to the agenda.

        The function f must except one argument, which is the the pointer to
        the workspace object on which the callback is executed.

        Parameters:

            f(callable): Python callable.

        """
        callback = c.CFUNCTYPE(None, c.c_void_p)(f)
        arts_api.agenda_insert_callback(self.ptr, callback)
        arts_api.callbacks += [callback]


    def execute(self, ws):
        """ Execute this agenda on the given workspace.
        Args:
            ws(Workspace): Workspace object on wich to execute the agenda.
        Raises:
            Exception: If execution of agenda on workspace fails.
        """
        with CoutCapture(ws):
            e = arts_api.execute_agenda(ws.ptr, self.ptr)
        if (e):
            raise Exception("Error during execution of Agenda:\n" + e.decode("utf8"))

    def __to_value_struct__(self):
        return {'ptr' : self.ptr}

    def __del__(self):
        """Destroys ARTS C API Agenda object associated with this Agenda object."""
        try:
            arts_api.destroy_agenda(self.ptr)
        except:
            pass

    @staticmethod
    def parse(name):
        """Parse controlfile and return agenda representing the agenda.

        Due to the way how ARTS works, agendas may not define WSVs with
        the same name. In this case, parsing of the agenda will fail.
        It is therefore not possible to parse an agenda twice which defines
        new workspace variables. In this case the user will have to keep track
        of the fist parsed agenda object.

        Args:
            name(str): Name of the control file. Is looked up recursively in the path
            specified by the ARTS_INCLUDE_PATH environmental variable.
        Raises:
            Exception: If parsing of the controlfile fails.
        """
        path = find_controlfile(name)
        ptr = arts_api.parse_agenda(path.encode())
        if not ptr:
            e = arts_api.get_error().decode("utf8")
            raise Exception("Error during parsing of controlfile " + str(path) + ":\n" + e)
        return Agenda(ptr)
