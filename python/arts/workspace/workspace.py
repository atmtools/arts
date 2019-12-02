"""
The workspace submodule.

Contains the Workspace which implements the main functionality of the ARTS interface.
Users should only have to use this class to interact with ARTS.

Attributes:
     imports(dict): Dictionary of parsed controlfiles. This is kept to ensure to avoid
                    crashing of the ARTS runtime, when a file is parsed for the second time.

"""
import ctypes as c
import logging
import numpy  as np

import ast
from   ast      import iter_child_nodes, parse, NodeVisitor, Call, Attribute, Name, \
                       Expression, Expr, FunctionDef, Starred, Module, expr
from   inspect  import getsource, getclosurevars
from contextlib import contextmanager
from copy       import copy
from functools  import wraps
import os

from arts.workspace.api       import arts_api, VariableValueStruct, \
                                            data_path_push, data_path_pop, \
                                            include_path_push, include_path_pop
from arts.workspace.methods   import WorkspaceMethod, workspace_methods
from arts.workspace.variables import WorkspaceVariable, group_names, group_ids, \
                                            workspace_variables
from arts.workspace.agendas   import Agenda
from arts.workspace import variables as V
from arts.workspace.output import CoutCapture
from arts.workspace.utility import unindent

imports = dict()


logger = logging.getLogger(__name__)

################################################################################
# ARTS Agenda Macro
################################################################################

class Include:
    """Simple helper class to handle INCLUDE statements in agenda definitions.

    Attributes:

        agenda: The included controlfile or agenda as
            arts.workspace.agenda.Agenda object.
    """
    def __init__(self, agenda):
        """ Create include from argument.

        Args:

            agenda (str, Agenda): Argument to the INCLUDE statement. This can
                either be a string or an Agenda object.
        """
        if type(agenda) == str:
            if not agenda in imports:
                self.agenda = Agenda.parse(agenda)
                imports[agenda] = self.agenda
            else:
                self.agenda = imports[agenda]
        elif type(agenda) == Agenda:
            self.agenda = agenda
        else:
            raise Exception("agenda argument must be either a controlfile"
            " name or a arts.workspace.agenda.Agenda object.")

def arts_agenda(func):
    """
    Parse python method as ARTS agenda

    This decorator can be used to define ARTS agendas using python function syntax.
    The function should have one arguments which is assumed to be a Workspace instance.
    All expressions inside the function must be calls to ARTS WSMs. The result is an
    Agenda object that can be used to copied into a named ARTS agenda

    Example:

    >>> @arts_agenda
    >>> def inversion_iterate_agenda(ws):
    >>>     ws.x2artsStandard()
    >>>     ws.atmfields_checkedCalc()
    >>>     ws.atmgeom_checkedCalc()
    >>>     ws.yCalc()
    >>>     ws.VectorAddVector(ws.yf, ws.y, ws.y_baseline)
    >>>     ws.jacobianAdjustAfterIteration()
    >>>
    >>> ws.Copy(ws.inversion_iterate_agenda, inversion_iterate_agenda)
    """

    source = getsource(func)
    source = unindent(source)
    ast = parse(source)

    func_ast = ast.body[0]
    if not type(func_ast) == FunctionDef:
        raise Exception("ARTS agenda definition can only decorate function definiitons.")

    args = func_ast.args.args

    try:
        arg_name = func_ast.args.args[0].arg
    except:
        raise Exception("Agenda definition needs workspace arguments.")

    ws = Workspace(0)

    context = copy(func.__globals__)
    context.update({arg_name : ws})
    # Add resolved non-local variables from closure.
    nls, _, _, _ = getclosurevars(func)
    context.update(nls)

    #
    # Helper functions
    #

    callback_body = []
    def callback_make_fun(body):
        """
        Helper function that creates a wrapper function around
        python code to be executed withing an ARTS agenda.
        """
        m = Module(body)

        def callback(ptr):
            try:
                context[arg_name].ptr = ptr
                eval(compile(m , "<unknown>", 'exec'), context)
            except Exception as e:
                logger.error(r"Exception in Python callback:\n", e)
            context[arg_name].ptr = None

        callback_body = []
        return callback

    def eval_argument(expr):
        """
        Evaluate argument of workspace method call.
        """
        if not hasattr(expr, "lineno"):
            setattr(expr, "lineno", 0)
        return eval(compile(Expression(expr), "<unknown>", 'eval'), context)

    # Create agenda
    a_ptr = arts_api.create_agenda(func.__name__.encode())
    agenda = Agenda(a_ptr)

    illegal_statement_exception = Exception(
        "Agenda definitions may only contain calls to WSMs of the"
        "workspace argument " + arg_name + " or INCLUDE statements.")

    #
    # Here the body of the function definition is traversed. Cases
    # that are treated specieal are INCLUDE statements and calls
    # of workspace methods. Remaining statements are accumulated
    # in callback_body and then added to the agenda as a single callback.
    #

    for e in func_ast.body:
        if not isinstance(e, Expr):
            callback_body += [e]
            continue
        else:
            call = e.value

        if not isinstance(call, Call):
            callback_body += [e]
            continue

        # Include statement
        if type(call.func) == Name:
            if not call.func.id == "INCLUDE":
                callback_body += [e]
            else:
                args = []
                for a in call.args:
                    args.append(eval_argument(a))
                    include = Include(*args)

                    if len(callback_body) > 0:
                        agenda.add_callback(callback_make_fun(callback_body))
                        callback_body = []

                    arts_api.agenda_append(agenda.ptr, include.agenda.ptr)
        else:
            att  = call.func.value
            if not att.id == arg_name:
                callback_body += [e]
                continue

            # Extract method name.
            name = call.func.attr

            # m is not a workspace method
            if not name in workspace_methods:
                callback_body += [e]
                continue

            # m is a workspace method.
            m  = workspace_methods[name]

            args = [ws, m]

            for a in call.args:
                # Handle starred expression
                if type(a) == Starred:
                    bs = eval_argument(a.value)
                    for b in bs:
                        args.append(b)
                    continue

                args.append(eval_argument(a))

            # Extract keyword arguments
            kwargs = dict()
            for k in call.keywords:
                kwargs[k.arg] = eval(
                    compile(Expression(k.value), "<unknown>", 'eval'),
                    context)

            # Add function to agenda
            if len(callback_body) > 0:
                agenda.add_callback(callback_make_fun(callback_body))
                callback_body = []

            agenda.add_method(*args, **kwargs)

    # Check if there's callback code left to add to the agenda.
    if len(callback_body) > 0:
        agenda.add_callback(callback_make_fun(callback_body))
        callback_body = []

    return agenda


################################################################################
# Workspace Method Wrapper Class
################################################################################
class WSMCall:
    """
    Wrapper class for workspace methods. This is necessary to be able to print
    the method doc as __repr__, which doesn't work for python function objects.

    Attributes:

        ws: The workspace object to which the method belongs.
        m:  The WorkspaceMethod object

    """
    def __init__(self, ws, m):
        self.ws = ws
        self.m  = m
        self.__doc__  = m.__doc__

    def __call__(self, *args, **kwargs):
        self.m.call(self.ws, *args, **kwargs)

    def __repr__(self):
        return repr(self.m)

################################################################################
# The Workspace Class
################################################################################
class Workspace:
    """
    The Workspace class represents an ongoing ARTS simulation. Each Workspace object
    holds its own ARTS workspace and can be used to execute ARTS workspace methods or
    access workspace variables.

    All workspace methods taken from workspace_methods in the methods module are added
    as attributed on creation and are thus available as class methods.

    Attributes:

        ptr(ctypes.c_void_p): object pointing to the ArtsWorkspace instance of the
            ARTS C API
        _vars(dict): Dictionary holding local variables that have been created
            interactively using the one of Create ARTS WSMs.


    """
    def __init__(self, verbosity=1, agenda_verbosity=0):
        """
        The init function just creates an instance of the ArtsWorkspace class of the
        C API and sets the ptr attributed to the returned handle.

        It also adds all workspace methods as attributes to the object.

        Parameters:
            verbosity (int): Verbosity level (0-3), 1 by default
            agenda_verbosity (int): Verbosity level for agendas (0-3),
                                    0 by default
        """

        self.__dict__["_vars"] = dict()
        self.ptr = arts_api.create_workspace(verbosity, agenda_verbosity)
        self.workspace_size = arts_api.get_number_of_variables()
        for name in workspace_methods:
            m = workspace_methods[name]
            setattr(self, m.name, WSMCall(self, m))
        self.__verbosity_init__()

    def __del__(self):
        """
        Cleans up the C API.
        """
        if not self.ptr is None:
            if not arts_api is None:
                arts_api.destroy_workspace(self.ptr)

    def __getstate__(self):
        return None

    def __setstate__(self):
        pass

    def __verbosity_init__(self):
        """
        Executes verbosityInit WSM directly through the ARTS api to suppress
        output.
        """
        wsm = workspace_methods["verbosityInit"]
        (m_id, args_out, args_in, ts) = wsm._parse_output_input_lists(self, [], {})
        arg_out_ptr = c.cast((c.c_long * len(args_out))(*args_out),
                             c.POINTER(c.c_long))
        arg_in_ptr  = c.cast((c.c_long * len(args_in))(*args_in),
                            c.POINTER(c.c_long))
        with CoutCapture(self, silent = True):
            e_ptr = arts_api.execute_workspace_method(self.ptr, m_id, len(args_out),
                                                      arg_out_ptr, len(args_in), arg_in_ptr)
        for t in ts[::-1]:
            t.erase()

    def create_variable(self, group, name):
        """
        Create a workspace variable.

        Args:

            group: The group name of the variable to create.

            name: The name of the variable to create. If None, the
            ARTS API will assign a unique name.

        """
        if not name is None:
            name = name.encode()

        group_id = group_ids[group]
        ws_id    = arts_api.add_variable(self.ptr, group_id, name)
        v        = arts_api.get_variable(ws_id)
        wsv      = WorkspaceVariable(ws_id,
                                     v.name.decode(),
                                     group_names[group_id],
                                     "User defined variable.",
                                     self)
        self._vars[wsv.name] = wsv
        return wsv

    def add_variable(self, var):
        """
        This will try to copy a given python variable to the ARTS workspace and
        return a WorkspaceVariable object representing this newly created
        variable.

        Types are natively supported by the C API are int, str, [str], [int], and
        numpy.ndarrays. These will be copied directly into the newly created WSV.

        In addition to that all arts types the can be stored to XML can
        be set to a WSV, but in this case the communication will happen through
        the file system (cf. WorkspaceVariable.from_arts).

        The user should not have to call this method explicitly, but instead it
        is used by the WorkspaceMethod call function to transfer python
        variable arguments to the ARTS workspace.

        Args:
            var: Python variable of type int, str, [str], [int] or np.ndarray
            which should be copied to the workspace.
        """
        if type(var) == WorkspaceVariable:
            return var

        # Create WSV in ARTS Workspace
        group = group_names[WorkspaceVariable.get_group_id(var)]
        wsv = self.create_variable(group, None)

        # Set WSV value using the ARTS C API
        s  = VariableValueStruct(var)
        if s.ptr:

            e = arts_api.set_variable_value(self.ptr, wsv.ws_id, wsv.group_id, s)
            if e:
                arts_api.erase_variable(self.ptr, wsv.ws_id, wsv.group_id)
                raise Exception("Setting of workspace variable through C API "
                                " failed with  the " + "following error:\n"
                                + e.decode("utf8"))
        # If the type is not supported by the C API try to write the type to XML
        # and read into ARTS workspace.
        else:
            try:
                wsv.from_arts(var)
            except:
                raise Exception("Could not add variable since + "
                                + str(type(var)) + " is neither supported by "
                                + "the C API nor arts XML IO.")
        self._vars[wsv.name] = wsv
        return wsv

    def __dir__(self):
        return {**self._vars, **workspace_variables, **self.__dict__}

    def __getattr__(self, name):
        """ Lookup the given variable in the local variables and the ARTS workspace.

        Args:
            name(str): Name of the attribute (variable)

        Raises:
            ValueError: If the variable is not found.
        """

        group_id = None
        if name in self._vars:
            var = self._vars[name]
            var.update()
            return var
        else:
            i = arts_api.lookup_workspace_variable(name.encode())
            if i < 0:
                raise AttributeError("No workspace variable " + str(name) + " found.")
            vs = arts_api.get_variable(i)
            group_id    = vs.group
            description = vs.description.decode("utf8")

        # Get its symbolic representation
        wsv = WorkspaceVariable(i, name, group_names[group_id], description, self)
        return wsv

    def __setattr__(self, name, value):
        """ Set workspace variable.

        This will lookup the workspace variable name and try to set it to value.

        Args:
            name(str):  Name of the attribute (variable)
            value(obj): The value to set the workspace variable to.

        Raises:
            ValueError: If the variable is not found or if value cannot uniquely converted to
            a value of a workspace variable.
        """
        try:
            v = self.__getattr__(name)
        except:
            self.__dict__[name] = value
            return None

        # Handle empty list or None values.
        if value is None or (isinstance(value, list) and not value):
            arts_api.set_variable_value(self.ptr, v.ws_id, v.group_id,
                                        VariableValueStruct.empty())
            return None

        if type(value) == Agenda:
            arts_api.set_variable_value(self.ptr, v.ws_id, v.group_id,
                                        VariableValueStruct(value))
            return None

        t = self.add_variable(value)

        if not t.group_id == v.group_id:
            raise Exception("Incompatible groups: Workspace variable " + name +
                            " of group " + group_names[v.group_id] + " and value " + str(value)
                            + " of group " + group_names[t.group_id] + ".")

        self.Copy(v, t)

        # Remove t only if it wasn't an existing WSV already before.
        if not type(value) == WorkspaceVariable:
            t.erase()

    def execute_agenda(self, agenda):
        """ Execute agenda on workspace.

        Args:

            agenda (arts.workspace.agenda.Agenda): Agenda object to execute.

        Raises:

            ValueError: If argument is not of type arts.workspace.agenda.Agenda
        """

        value_error = ValueError("Argument must be of type agenda.")
        if not type(agenda) is Agenda:
            raise value_error

        include_path_push(os.getcwd())
        data_path_push(os.getcwd())

        agenda.execute(self)

        include_path_pop()
        data_path_pop()

    def execute_controlfile(self, name):
        """ Execute controlfile or agenda on workspace.

        This method looks recursively for a controlfile with the given name in the current
        directory and the arts include path. If such a file has been found it will be parsed
        and executed on the workspace.

        Args:

            name(str): Name of the controlfile

        Raises:

            Exception: If parsing of the controlfile fails.

        Returns:

            The controlfile as parsed Agenda object.

        """

        if not name in imports:
            agenda = Agenda.parse(name)
            imports[name] = agenda
        else:
            agenda = imports[name]

        self.execute_agenda(agenda)

        return agenda
