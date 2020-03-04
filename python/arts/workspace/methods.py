""" The methods submodule.

This module exposes all ARTS workspace methods represented by WorkspaceMethod object.

The methods are loaded dynamically when the module is imported, which ensures that they
up to date with the current ARTS build.

Attributes:

     workspace_methods(dict): Dictionary containing all ARTS workspace methods.

"""

import ast
import re
import ctypes as c
import numpy as np


from arts.workspace.api       import arts_api, nodef
from arts.workspace.variables import WorkspaceVariable, group_ids, group_names
from arts.workspace import variables, workspace
from arts.workspace.output import CoutCapture

class WorkspaceMethod:
    """
    The WorkspaceMethod class represents ARTS workspace methods. Each workspace method
    provided a call function that forwards the call to the function towards the ARTS C API.

    Attributes:
        m_ids([int]):       Indices of supergeneric overloads of this WSM
        name(str):          The name of the method as defined in methods.cc
        description(str):   The documentation of the method as defined in methods.cc
        outs([int]):        Indices of the output variables of the method.
        n_out(int):         The number of output arguments.
        n_g_out(int):       The number of generic outputs.
        g_out([str]):       The names of the generic output arguments.
        g_out_types([dict]): List of dicts associating the name of a generic output
                             with its types for each of the supergeneric overloads.
        n_in(int):          The number of input arguments.
        ins([int]):         The indices of the input arguments of the WSM.
        n_g_in(int):        The number of generic input arguments.
        g_in_types([dict]): List of dicts associating the name of a generic input to the
                            expected type for each supergeneric overload of the method.
        g_in_default(dict)  Dictionary containing the default values for each generic parameter.
        g_in([str]):        List of names of the generic input arguments.
    """

    # Regular expression to that matches <Group>Create WSMs
    create_regexp = re.compile(r"^(\w*)Create$")

    def __init__(self, m_id, name, description, outs, g_out_types, ins, g_in_types):
        """Create a WorkspaceMethod object from a given id, name, description and types of
        generic input and output arguments.

        Args:
        m_id(int):          The index identifying the method in the C API
        name(str):          The name of the method
        description(str):   Method documentation
        outs([int]):        List of indices of the (non-generic) output arguments of the method.
        g_out_types([int]): Group ids of the generic output types
        outs([int]):        List of indices of the (non-generic) input arguments of the method.
        g_in_types([int]):  Group ids of the generic input types
        """
        self.m_ids       = [m_id]
        self.name        = name
        self.description = description

        # Output
        self.outs  = outs
        self.n_out = len(outs)

        # Generic Output
        self.n_g_out     = len(g_out_types)
        self.g_out_types = [WorkspaceMethod.get_output_dict(m_id, g_out_types)]
        self.g_out       = [k for k in self.g_out_types[0]]

        # Input
        self.ins  = ins
        self.n_in = len(ins)

        # Generic Input
        self.n_g_in       = len(g_in_types)
        self.g_in_types   = [WorkspaceMethod.get_input_dict(m_id, g_in_types)]
        self.g_in_default = WorkspaceMethod.get_default_input_dict(m_id, g_in_types)
        self.g_in         = [k for k in self.g_in_types[0]]

        self.is_create = False
        if (WorkspaceMethod.create_regexp.match(name)):
                self.is_create = True

    def __repr__(self):
        return "ARTS Workspace Method: " + self.name

    def help(self):
        """
        Print ARTS documentation for this method.
        """
        print(arts_api.method_print_doc(self.m_ids[0]).decode("utf8"))
        return None

    def add_overload(self, m_ids, g_in_types, g_out_types):
        """ Add one or more overloads to a workspace method.

        Use this function to add a supergeneric overload to a WorkspaceMethod object
        so that is will be considered in overload resolution when call(...) is called.

        TODO: Simplify this to take a WorkspaceMethod object.

        Args:
            m_ids([int]): The method ids of the supergeneric ARTS WSMs which should be added
                          to the list of overloads.
            g_in_types ([dict]): List of dicts containing the mappings between argument
                                 names and indices of expected groups
            g_out_types ([dict]): List of dicts containing the mappings between argument
                                  names and indices of expected groups
        """
        self.m_ids += m_ids
        self.g_in_types += g_in_types
        self.g_out_types += g_out_types

    @staticmethod
    def get_input_dict(m_id, in_types):
        """Get mapping of names of generic input variables to indices of groups.

        Args:
            m_id(int):       Index of the method.
            in_types([int]): List of group indices of the generic input arguments.
        Return:
            dict: The mapping.
        """
        res = dict()
        for i,t in enumerate(in_types):
            res[arts_api.get_method_g_in(m_id, i).decode("utf8")] = t
        return res

    @staticmethod
    def get_default_input_dict(m_id, g_in):
        """Get dict mapping names of generic input arguments to default values.

        Is None if no default value for a given generic input is given.

        Args:
            m_id(int):   Index of the method.
            g_in([str]): Names of the generic input arguments.
        Return:
            dict: The mapping.
        """
        res = dict()
        for i,t in enumerate(g_in):
            k = arts_api.get_method_g_in(m_id, i).decode("utf8")
            d = arts_api.get_method_g_in_default(m_id, i).decode("utf8").strip()
            if d == nodef:
                d = None
            elif d == "Inf":
                res[k] = np.float64("inf")
            else:
                try:
                    d = WorkspaceVariable.convert(group_names[t], ast.literal_eval(d))
                    res[k] = d
                except:
                    res[k] = d
        return res

    @staticmethod
    def get_output_dict(m_id, out_types):
        """Get mapping of names of generic output variables to indices of groups.

        Args:
            m_id(int):       Index of the method.
            in_types([int]): List of group indices of the generic output arguments.
        Return:
            dict: The mapping.
        """
        res = dict()
        for i,t in enumerate(out_types):
            res[arts_api.get_method_g_out(m_id, i).decode("utf8")] = t
        return res

    def _parse_output_input_lists(self, ws, args, kwargs):

        ins  = [i for i in self.ins if not i in self.outs]
        outs = self.outs[:]
        temps = []

        # Can call function using only positional or named arguments
        if len(args) and len(kwargs.keys()):
            raise SyntaxError("ARTS WSMs can only be called using either " +
                              " positional or named arguments.")

        # Use dict to store named arguments. Generic in- and outputs
        # will be added to this.
        named_args = dict(**kwargs)

        #
        # Parse positional arguments
        #

        for j in range(len(args)):
            # Parse output arguments
            if j < self.n_out:
                if not type(args[j]) == WorkspaceVariable:
                    out_name = WorkspaceVariable.get_variable_name(outs[j])
                    raise TypeError("Positional argument for output {} is not of type "
                                    "WorkspaceVariable.".format(out_name))

                outs[j] = args[j].ws_id
            # Parse generic output arguments and to list of named
            # args.
            elif j < self.n_out + self.n_g_out:
                k = j - self.n_out
                name = self.g_out[k]
                named_args[name] = args[j]

            # Parse input arguments
            elif j < self.n_out + self.n_g_out + len(ins):
                k = j - self.n_g_out - self.n_out
                if type(args[j]) == WorkspaceVariable:
                    ins[k] = args[j].ws_id
                else:
                    temps.append(ws.add_variable(args[j]))
                    ins[k] = temps[-1].ws_id
            # Parse generic input arguments
            elif j < self.n_out + self.n_g_out + len(ins) + self.n_g_in:
                k = j - self.n_g_out - self.n_out - len(ins)
                name = self.g_in[k]
                named_args[name] = args[j]
            else:
                raise Exception(" {} positional arguments given, but this WSM takes at "
                                "most {}.".format(len(args), j))

        # Parse named arguments
        ins_names = [WorkspaceVariable.get_variable_name(i) for i in ins]
        outs_names = [WorkspaceVariable.get_variable_name(i) for i in outs]

        # Raise exception if named argument does not match method arguments.
        for k in kwargs:
            if not (k in ins_names or k in outs_names or
                    k in self.g_in or k in self.g_out):
                raise Exception("The provided named argument '{0}' does not match "
                                " any of the arguments of WSM {1}."
                                .format(k, self.name))

            if k in ins_names:
                i = ins_names.index(k)
                arg = kwargs[k]
                if type(arg) == WorkspaceVariable:
                    ins[i] = arg.ws_id
                else:
                    temps.append(ws.add_variable(arg))
                    ins[i] = temps[-1].ws_id

            if k in outs_names:
                i = outs_names.index(k)
                arg = kwargs[k]
                if type(arg) == WorkspaceVariable:
                    outs[i] = arg.ws_id
                else:
                    raise Exception(("Output argument {} must be a workspace "
                                     + "variable.").format(k))

        # Check output argument names
        g_output_args = dict()
        for k in self.g_out:
            if not k in named_args:
                raise Exception("WSM " + self.name + " needs generic output " + k)
            else:
                g_output_args[k] = named_args[k]

        # Check input argument names
        g_input_args = dict()
        for k in self.g_in:
            if not k in named_args:
                if k in self.g_in_default:
                    g_input_args[k] = self.g_in_default[k]
                else:
                    raise Exception("WSM " + self.name + " needs generic input " + k)
            else:
                g_input_args[k] = named_args[k]

        m_id = self.m_ids[0]
        sg_index = 0

        if (len(self.m_ids) > 1):

            #
            # First, try resolving solely based on WSV input since we can
            # easily identify their type.
            #

            candidates = []
            for i, m_id in enumerate(self.m_ids):
                g_in_types = self.g_in_types[i]
                g_out_types = self.g_out_types[i]

                convertable = True

                for k in g_in_types:
                    in_type = group_names[g_in_types[k]]
                    in_arg = g_input_args[k]
                    if isinstance(in_arg, WorkspaceVariable):
                        if not in_arg.group == in_type:
                            convertable = False
                            break

                for k in g_out_types:
                    out_type = group_names[g_out_types[k]]
                    out_arg = g_output_args[k]
                    if isinstance(out_arg, WorkspaceVariable):
                        if not out_arg.group == out_type:
                            convertable = False
                            break

                if convertable:
                    candidates.append((m_id, i))

            if len(candidates) == 1:
                m_id, sg_index = candidates[0]
            else:
                # Resolve overload (if necessary).
                g_out_types = dict([(k,WorkspaceVariable.get_group_id(g_output_args[k]))
                                    for k in self.g_out])
                g_in_types  = dict([(k,WorkspaceVariable.get_group_id(g_input_args[k]))
                                    for k in self.g_in])

                out_indices = [i for i,ts in enumerate(self.g_out_types) if ts == g_out_types]
                in_indices  = [i for i,ts in enumerate(self.g_in_types) if ts == g_in_types]
                sg_indices  = set(out_indices) & set(in_indices)

                if len(sg_indices) > 1:
                    raise Exception("Could not uniquely resolve super-generic overload.")

                if len(sg_indices) == 0:
                    raise Exception("Could not find super-generic overload matching"
                                    + " the given groups.")

                sg_index = sg_indices.pop()
                m_id = self.m_ids[sg_index]

        # Combine input and output arguments into lists.
        arts_args_out = []
        for out in outs:
            arts_args_out.append(out)

        for name in self.g_out:
            arg = g_output_args[name]
            if not type(arg) == WorkspaceVariable:
                raise ValueError("Generic Output " + name + " must be an ARTS WSV.")
            group_id = arg.group_id
            expected = self.g_out_types[sg_index][name]
            if not group_id == expected:
                raise Exception("Generic output " + name + " expected to be of type "
                                + group_names[expected])
            arts_args_out.append(arg.ws_id)

        arts_args_in = []
        for i in ins:
            if not i in outs:
                arts_args_in.append(i)

        for name in self.g_in:
            arg = g_input_args[name]
            if type(arg) == WorkspaceVariable:
                arts_args_in.append(arg.ws_id)
            else:
                gid = self.g_in_types[sg_index][name]
                arg_converted = WorkspaceVariable.convert(group_names[gid], arg)
                if arg_converted is None:
                    raise Exception("Could not convert input {} to expected group {}"
                                    .format(arg, group_names[gid]))
                temps.append(ws.add_variable(arg_converted, gid))
                arts_args_in.append(temps[-1].ws_id)
        return (m_id, arts_args_out, arts_args_in, temps)


    def create(self, ws, name = None):
        """
        Call to <Group>Create WSMs are handled differently. This method simply
        determines the group type from the function name and then add a variable of
        this type to the workspace ws. A handle of this variable is then added to
        as attribute to the arts.workspace.variables module.

        Args:
            ws(Workspace): Workspace object to add the variable to
            name(str):     Name of the variable to add to the workspace
        """
        group = WorkspaceMethod.create_regexp.match(self.name).group(1)
        group_id = group_ids[group]

        if not name:
            name = "__anonymous_" + str(len(ws._vars))
            ws_id = arts_api.add_variable(ws.ptr, group_id, name.encode())
        else:
            # Is there a WSM with that name?
            if name in workspace_methods.keys():
                raise Exception("A WSM with the name " + name + " already exists.")

            # Is there a WSV with that name?
            ws_id = arts_api.lookup_workspace_variable(name.encode())
            # Is yes, check that it is of the same group?
            if not ws_id == -1:
                v = arts_api.get_variable(ws_id)
                if not v.group == group_id:
                    raise Exception("A WSV with the name " + name + " but of goup "
                                    + group_names[v.group] + " already exists.")
            # Otherwise we add the variable.
            else:
                ws_id = arts_api.add_variable(ws.ptr, group_id, name.encode())

        wsv = WorkspaceVariable(ws_id, name, group, "User defined variable.", ws)
        setattr(variables, name, wsv)
        ws._vars[name] = wsv
        return wsv

    def call(*args, **kwargs):
        """ Execute workspace method.

        This method will execute the workspace method (args[0]) on the workspace object (args[1])
        interpreting the remaining arguments in `*args` and `**kwargs` as arguments.

        Positional arguments in `*args` are interpreted in order with output arguments coming
        first.

        Keyword arguments in kwargs are interpreted according to the name of the generic
        parameters of the ARTS WSM.

        Args:
        args(list): Positional arguments with the first argument being the WorkspaceMethod
        instance, i.e. self = args[0], the second the Workspace object (args[1]). The
        remaining arguments are interpreted as generic arguments to the ARTS WSM.
        kargs(dict): Keyword args are interpreted as named generic arguments to the ARTS WSM
        according to its definition in methods.cc.
        """

        self = args[0]

        if self.is_create:
            return self.create(*args[1:])

        ws   = args[1]

        (m_id, arts_args_out, arts_args_in, temps) = self._parse_output_input_lists(ws,
                                                                                    args[2:],
                                                                                    kwargs)

        # Execute WSM and check for errors.
        arg_out_ptr = c.cast((c.c_long * len(arts_args_out))(*arts_args_out), c.POINTER(c.c_long))
        arg_in_ptr = c.cast((c.c_long * len(arts_args_in))(*arts_args_in), c.POINTER(c.c_long))

        with CoutCapture(ws):
            e_ptr = arts_api.execute_workspace_method(ws.ptr, m_id,
                                                      len(arts_args_out),
                                                      arg_out_ptr,
                                                      len(arts_args_in),
                                                      arg_in_ptr)
        if (e_ptr):
            raise Exception("Call to ARTS WSM " + self.name + " failed with error: "
                           + e_ptr.decode("utf8").format())

        # Remove temporaries from workspace (in reverse order).
        for t in temps[::-1]:
            t.erase()


    def describe(self):
        """ Print WSM documentation. """
        print(self.description.format())

def iter_raw():
    """ Iterator returning a WorkspaceMethod object for each available ARTS WSM.

    This iterator returns super-generically overloaded methods several times.

    Yields:
        WorkspaceMethod: The next ARTS Workspace method as defined in methods.cc in
                         increasing order.
    """
    for i in range(arts_api.get_number_of_methods()):
        m            = arts_api.get_method(i)
        name         = m.name.decode("utf8")
        description  = m.description.decode("utf8")
        outs         = [m.outs[i] for i in range(m.n_out)]
        g_out_types  = [m.g_out_types[i] for i in range(m.n_g_out)]
        ins          = [m.ins[i] for i in range(m.n_in)]
        g_in_types   = [m.g_in_types[i] for i in range(m.n_g_in)]
        yield WorkspaceMethod(m.id, name, description, outs, g_out_types, ins, g_in_types)

def iter():
    """ Iterator returning a WorkspaceMethod object for each available ARTS WSM.

    This iterator returns overloaded Workspace methods, i.e. super-generically overloaded
    WSM are not returned multiple times.

    Yields:
        WorkspaceMethod: The next ARTS Workspace method as defined in methods.cc in
                         increasing order.
    """
    for k,m in workspace_methods:
        yield m

workspace_methods = dict()
for m in iter_raw():
    if m.name in workspace_methods:
        workspace_methods[m.name].add_overload(m.m_ids, m.g_in_types, m.g_out_types)
    else:
        workspace_methods[m.name] = m
