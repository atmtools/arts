import os
import sys
from   ast      import parse, Call, Name, Expression, Expr, FunctionDef, \
                       Starred, Module, Str
from   inspect  import getsource, getclosurevars
from   copy     import copy

import pyarts.pyarts_cpp as cxx
from pyarts.workspace.utility import unindent as unindent


# Set the default basename of Arts
try:
    filename = sys.modules["__main__"].__file__
    basename, _ = os.path.splitext(os.path.basename(filename))
    cxx.parameters.out_basename = basename
except:
    pass


InternalWorkspace = getattr(cxx, "Pyarts::Workspace")
Agenda = cxx.Agenda


class DelayedAgenda:
    """ Helper class to delay the parsing of an Agenda until a workspace exist
    """
    
    def __init__(self, *args):
        self.args = [[*args]]
    def append_agenda_methods(self, other):
        self.args.extend(other.args)
    def __call__(self, ws):
        a = cxx.Agenda()
        for args in self.args:
            a.append_agenda_methods(continue_parser_function(ws, *args))
        a.name = "<unknown>"
        return a


def Include(ws, path):
    """ Parse and execute the .arts file at path onto the current workspace ws
    
    The Arts parser is invoked on the file at path.  The methods and commands
    of this file are executed
    """
    if isinstance(path, Agenda):
        path.execute(ws)
    else:
        Agenda(ws, path).execute(ws)

def arts_agenda(func=None, *, ws=None, allow_callbacks=False, set_agenda=False):
    """
    Decorator to parse a Python method as ARTS agenda

    This decorator can be used to define ARTS agendas using python function
    syntax. The function should have one arguments which is assumed to be a
    Workspace instance. All expressions inside the function must be calls to
    ARTS WSMs. The function definition results in an Agenda object that can
    be copied into an ARTS agenda.

    Example:

    >>> @arts_agenda(ws)
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

    >>> @arts_agenda(ws, allow_callbacks=True)
    >>> def python_agenda(ws):
    >>>     print("Python says 'hi'.)
                  
    Note that using allow_callbacks=True breaks the Agenda input-output
    control.  It is therefore considered undefined behavior if you manipulate
    workspace variables that are neither in- nor output of the Agenda using
    callbacks.  Do this at your own risk.
                  
    A special INCLUDE(path) directive can be part of the function definitions
    to use the Arts parser of .arts files to be invoked on the file.
    All methods and invokations that are part of the .arts file are
    appended in place to the agenda
    
    If set_agenda is set True, the workspace agenda is set to the function name
    (allowing initiating automatic checking if the type is named as a defined
    workspace agenda)
    """
    
    def agenda_decorator(func):
        return parse_function(func, ws, allow_callbacks=allow_callbacks, set_agenda=set_agenda)
    
    if func is None:
        return agenda_decorator
    else:
        return parse_function(func, ws, False, False)

def parse_function(func, arts, allow_callbacks, set_agenda):
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
    
    context = copy(func.__globals__)
    nls, _, _, _ = getclosurevars(func)
    context.update(nls)
    
    if arts is None:
        return DelayedAgenda(context, ast, allow_callbacks, set_agenda)
    return continue_parser_function(arts, context, ast, allow_callbacks, set_agenda)


def continue_parser_function(arts, context, ast, allow_callbacks, set_agenda):
    assert isinstance(arts, InternalWorkspace), f"Expects Workspace, got {type(arts)}"

    func_ast = ast.body[0]
    if not isinstance(func_ast, FunctionDef):
        raise Exception("ARTS agenda definition can only decorate function definitions.")
    
    args = func_ast.args.args

    try:
        arg_name = func_ast.args.args[0].arg
    except:
        raise Exception("Agenda definition needs workspace arguments.")
    
    context.update({arg_name : arts})
    
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
        
        def callback(ws):
            # FIXME: Is mutex required here?
            context[arg_name].swap(ws)  # FIXME:  This is required, of course
            eval(compile(m , "<unknown>", 'exec'), context)
            context[arg_name].swap(ws)  # FIXME: But is this required or thrown away?
        return callback

    def eval_argument(expr):
        """
        Evaluate argument of workspace method call.
        """
        if not hasattr(expr, "lineno"):
            setattr(expr, "lineno", 0)
        return eval(compile(Expression(expr), "<unknown>", 'eval'), context)

    illegal_statement_exception = Exception(
        "Pure ARTS agenda definitions may only contain calls to WSMs of"
        " the workspace argument '{arg_name}' or INCLUDE statements."
        " If you want to allow Python callbacks you need to use"
        " the '@arts_agenda' decorator with the 'allow_callbacks'"
        " keyword argument set to 'True'.")

    workspace_methods = [str(x.name) for x in cxx.get_md_data()]
    agenda = Agenda()
    
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
                    include_agenda =  Agenda(arts, *args)

                    if len(callback_body) > 0:
                        agenda.add_callback_method(arts, callback_make_fun(callback_body))
                        callback_body = []

                    agenda.append_agenda_methods(include_agenda)
        else:
            att = call.func.value
            if not hasattr(att, 'id') or not att.id == arg_name:
                callback_body += [e]
                continue

            # Extract method name.
            name = call.func.attr
            # m is not a workspace method
            if name not in workspace_methods:
                if allow_callbacks:
                    callback_body += [e]
                    continue
                else:
                    raise ValueError(
                        f"{name} is not a know ARTS WSM."
                    )

            args = []
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
                agenda.add_callback_method(arts, callback_make_fun(callback_body))
                callback_body = []
            
            agenda.add_workspace_method(arts, name, *args, **kwargs)

    # Check if there's callback code left to add to the agenda.
    if len(callback_body) > 0:
        agenda.add_callback_method(arts, callback_make_fun(callback_body))
        callback_body = []
    
    agenda.name = func_ast.name
    if set_agenda: setattr(arts, func_ast.name, agenda)
    return agenda


_group_types = [type(eval(f"cxx.{x}()")) for x in list(cxx.get_wsv_group_names())]


class Workspace(InternalWorkspace):
    def __getattr__(self, attr):
        if self._hasattr_check_(attr): return self._getattr_unchecked_(attr)
        if attr == "__class__": return InternalWorkspace
        raise AttributeError(f"'Workspace' object has no attribute '{attr}'")
    
    def __setattr__(self, attr, value):
        if self._hasattr_check_(attr):
            if isinstance(value, DelayedAgenda): value = value(self)
            self._getattr_unchecked_(attr).initialize_if_not()
            self._getattr_unchecked_(attr).value = type(self._getattr_unchecked_(attr).value)(value)
        else:
            if type(value) in _group_types:
                self._setattr_unchecked_(attr, value)
            elif isinstance(value, DelayedAgenda):
                self._setattr_unchecked_(attr, value(self))
            else:
                raise AttributeError(f"'Workspace' object has no attribute '{attr}'")
    
