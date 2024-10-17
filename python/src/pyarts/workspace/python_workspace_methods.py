from .workspace import python_workspace_methods
from ast import parse, FunctionDef, Return, Name, Tuple, List, ClassDef
from pyarts.workspace.utility import unindent as unindent
from inspect import getsource
import  pyarts.arts as cxx


class InnerPythonWorkspaceMethod:
    def __init__(self, input, output, func):

        for out in output:
            assert out in cxx.globals.workspace_variables(), f"Output variable {out} is not a workspace variable"
        for inp in input:
            assert inp in cxx.globals.workspace_variables(), f"Input variable {inp} is not a workspace variable, consider making it a keyword-only argument"

        self.input = input
        self.output = output
        self.func = func
      
    def __repr__(self):
        return f"InnerPythonWorkspaceMethod(in:, {self.input}, out: {self.output})"
    __str__ = __repr__

    def __call__(self, ws, *args, **kwargs):
        inner_kwargs = {}
        for n, a in zip(self.input, args):
            inner_kwargs[n] = a

        inner_kwargs.update(kwargs)

        for n in self.input:
            if n not in inner_kwargs:
                assert ws.has(n), f"Workspace variable {n} is not initialized"
                inner_kwargs[n] = getattr(ws, n)

        retval = self.func(**inner_kwargs)

        if len(self.output) == 1:
            if self.output[0] in inner_kwargs:
                if hasattr(inner_kwargs[self.output[0]], "_set_value"):
                    inner_kwargs[self.output[0]]._set_value(retval)
                else:
                    raise RuntimeError(f"Input/output argument for '{self.output[0]}' does not appear to be a workspace group, is type {type(inner_kwargs[self.output[0]])}")
            else:
                setattr(ws, self.output[0], retval)
        else:
            for n, v in zip(self.output, retval):
                if n in inner_kwargs:
                    if hasattr(inner_kwargs[n], "_set_value"):
                        inner_kwargs[n]._set_value(v)
                    else:
                      raise RuntimeError(f"Input/output argument for '{n}' does not appear to be a workspace group, is type {type(v)}")
                else:
                    setattr(ws, n, v)


def _callback_operator_return(expr):
    value = expr.value

    if isinstance(value, Name):
        return [value.id]
    elif isinstance(value, Tuple) or isinstance(value, List):
        return [x.id if isinstance(x, Name) else _errvar for x in value.elts]
    return None


_errvar = "Unknown Value"


def _find_return(body):
    bad = False
    for expr in body:
        if isinstance(expr, Return):
            x = _callback_operator_return(expr)
            if x is not None:
                return x
            else:
                bad = True
        elif isinstance(expr, FunctionDef) or isinstance(expr, ClassDef):
            continue

        if hasattr(expr, "body"):
            t = _find_return(expr.body)
            if t is not None:
                return t
        if hasattr(expr, "orelse"):
            t = _find_return(expr.orelse)
            if t is not None:
                return t

    if bad:
        return [_errvar]

    return []


def _callback_operator_function(funcdef):
    assert isinstance(funcdef, FunctionDef), "Must be a function"
    return [x.arg for x in funcdef.args.args]


def _python_workspace_method(func):
    """Internal source code parser"""
    srccod = getsource(func)
    srccod = unindent(srccod)
    srcast = parse(srccod)

    assert len(srcast.body) == 1

    code_body = srcast.body[0]

    args = _callback_operator_function(code_body)
    retvals = _find_return(code_body.body)
    ret = InnerPythonWorkspaceMethod(args, retvals, func)

    assert func.__name__ not in python_workspace_methods, "Function already exists"
    python_workspace_methods[func.__name__] = ret

    return ret


def python_workspace_method(func=None):
    """
    Creates a python workspace method
    """

    def parser(fn):
        return _python_workspace_method(fn)

    if func is None:
        return parser
    else:
        return _python_workspace_method(func)
