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
import sys
import numpy  as np
import weakref

import ast
from   ast      import iter_child_nodes, parse, NodeVisitor, Call, Attribute, Name, \
                       Expression, Expr, FunctionDef, Starred, Module, expr, Str
from   inspect  import getsource, getclosurevars
from contextlib import contextmanager
from copy       import copy
from functools  import wraps
import os

from pyarts.workspace.api import (arts_api, VariableValueStruct, data_path_push,
                                data_path_pop, include_path_push,
                                include_path_pop, is_empty)
from pyarts.workspace.methods   import WorkspaceMethod, workspace_methods
from pyarts.workspace.variables import (WorkspaceVariable, group_names, group_ids,
                                      workspace_variables)
from pyarts.workspace.agendas   import Agenda
from pyarts.workspace import variables as V
from pyarts.workspace.output import CoutCapture
from pyarts.workspace.utility import unindent

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


def arts_agenda(func=None, *, allow_callbacks=False):
    """
    Decorator to parse a Python method as ARTS agenda.

    This decorator can be used to define ARTS agendas using python function
    syntax. The function should have one arguments which is assumed to be a
    Workspace instance. All expressions inside the function must be calls to
    ARTS WSMs. The function definition results in an Agenda object that can
    be copied into an ARTS agenda.

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

    When the decorator is used with the 'allow_callbacks' keyword argument
    set to True, arbitrary Python code can be executed within the callback.
    Note, however, that ARTS ignores exceptions occurring within the
    callback, so care must be taken that potentially silenced errors
    don't interfere with simulation results.

    Example:

    >>> @arts_agenda(allow_callbacks=True)
    >>> def python_agenda(ws):
    >>>     print("Python says 'hi'.)
    """
    def agenda_decorator(func):
        return parse_function(func, allow_callbacks=allow_callbacks)

    if func is None:
        return agenda_decorator
    else:
        return parse_function(func, False)


def parse_function(func, allow_callbacks):
    """
    Parse python method as ARTS agenda

    Args:
        func: The function object to parse.
        allow_callbacks: Whether to allow callbacks in the agenda.

    Return:
        An 'Agenda' object containing the code in the given function.
    """

    source = getsource(func)
    source = unindent(source)
    ast = parse(source)
    print(func, "\n\n\n", type(source), "\n\n\n", ast)

    func_ast = ast.body[0]
    if not type(func_ast) == FunctionDef:
        raise Exception("ARTS agenda definition can only decorate function definitions.")

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
        if sys.version_info >= (3, 8):
            # https://bugs.python.org/issue35894#msg334808
            m = Module(body, [])
        else:
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
    print(func.__name__.encode())
    agenda = Agenda(a_ptr)

    illegal_statement_exception = Exception(
        "Pure ARTS agenda definitions may only contain calls to WSMs of"
        " the workspace argument '{arg_name}' or INCLUDE statements."
        " If you want to allow Python callbacks you need to use"
        " the '@arts_agenda' decorator with the 'allow_callbacks'"
        " keyword argument set to 'True'."
    )

    #
    # Here the body of the function definition is traversed. Cases
    # that are treated specieal are INCLUDE statements and calls
    # of workspace methods. Remaining statements are accumulated
    # in callback_body and then added to the agenda as a single callback.
    #

    for e in func_ast.body:
        if not isinstance(e, Expr):
            if allow_callbacks:
                callback_body += [e]
                continue
            else:
                raise illegal_statement_exception
        else:
            call = e.value

        if not isinstance(call, Call):
            if isinstance(call, Str):
                continue
            elif allow_callbacks:
                callback_body += [e]
                continue
            else:
                raise illegal_statement_exception

        # Include statement
        if type(call.func) == Name:
            if call.func.id != "INCLUDE":
                if allow_callbacks:
                    callback_body += [e]
                else:
                    raise illegal_statement_exception
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
            att = call.func.value
            if not att.id == arg_name:
                callback_body += [e]
                continue

            # Extract method name.
            name = call.func.attr
            print(workspace_methods)
            # m is not a workspace method
            if name not in workspace_methods:
                if allow_callbacks:
                    callback_body += [e]
                    continue
                else:
                    raise ValueError(
                        f"{name} is not a know ARTS WSM."
                    )

            # m is a workspace method.
            m = workspace_methods[name]

            args = [ws, m]
            kwargs = dict()

            for a in call.args:
                # Handle starred expression
                if type(a) == Starred:
                    bs = eval_argument(a.value)
                    for b in bs:
                        args.append(b)
                else:
                    args.append(eval_argument(a))

            # Extract keyword arguments
            for k in call.keywords:
                if k.arg is None:
                    d = eval(compile(Expression(k.value), "<unknown>", 'eval'),
                             context)
                    kwargs.update(d)
                else:
                    kwargs[k.arg] = eval(compile(Expression(k.value),
                                                 "<unknown>", 'eval'),
                                         context)

            # Add function to agenda
            if len(callback_body) > 0:
                agenda.add_callback(callback_make_fun(callback_body))
                callback_body = []
            print(args, '\n', kwargs)
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
        self._ws = weakref.ref(ws)
        self.m  = m
        self.__doc__  = m.__doc__

    @property
    def ws(self):
        return self._ws()

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
    def __init__(self, verbosity=0, agenda_verbosity=0):
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
                self.ptr = None

    def __getstate__(self):
        raise Exception("ARTS workspaces cannot be pickled.")

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

    def add_variable(self, var, group = None):
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
        if group is None:
            group = WorkspaceVariable.get_group_id(var)

        group = group_names[group]
        wsv = self.create_variable(group, None)

        # Set WSV value using the ARTS C API
        self.set_variable(wsv, var)
        self._vars[wsv.name] = wsv
        return wsv

    def set_variable(self, wsv, value):
        """
        This will set a WSV to the given value.

        Natively supported types, i.e. any of int, str, [str], [int], numpy.ndarrays,
        and scipy.sparse, will be copied directly into the newly created WSV.

        In addition to that all arts types the can be stored to XML can
        be set to a WSV, but in this case the communication will happen through
        the file system (cf. :code:`WorkspaceVariable.from_arts`).

        Args:
            wsv: The :class:`WorkspaceVariable` to set.
            value: The Python object representing the value to set :code:`wsv` to.
        """
        if is_empty(value):
            err = arts_api.set_variable_value(self.ptr, wsv.ws_id, wsv.group_id,
                                              VariableValueStruct.empty())
            if not err is None:
                msg = ("The following error occurred when trying to set the"
                       " WSV {}: {}".format(wsv.name, err.decode()))
                raise Exception(msg)
            return None

        group = group_names[WorkspaceVariable.get_group_id(value)]
        if group != wsv.group:
            try:
                converted = WorkspaceVariable.convert(wsv.group, value)
            except:
                converted = None
            if converted is None:
                raise Exception("Cannot set workspace variable of type {} to "
                                " value  '{}'.".format(wsv.group, value))
            value = converted

        s = VariableValueStruct(value)
        if s.ptr:
            err = arts_api.set_variable_value(self.ptr, wsv.ws_id, wsv.group_id, s)
            if not err is None:
                msg = ("The following error occurred when trying to set the"
                       " WSV {}: {}".format(wsv.name, err.decode()))
                raise Exception(msg)
        # If the type is not supported by the C API try to write the type to XML
        # and read into ARTS workspace.
        else:
            try:
                wsv.from_arts(value)
            except:
                raise Exception("Could not set variable since + "
                                + str(type(value)) + " is neither supported by "
                                + "the C API nor arts XML IO.")

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

        self.set_variable(v, value)

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
