from ast import parse, FunctionDef, Return, Name, Tuple, List, ClassDef
from inspect import getsource
import pyarts.arts as cxx
from pyarts.workspace import Workspace
from pyarts.workspace.utility import unindent as unindent


_errvar = "Unknown Value"


_group_types = [
    eval(f"cxx.{x.name}") for x in list(cxx.globals.get_wsv_groups())
]


def safe_set(ws, attr, val):
    if ws._hasattr_check_(attr):
        ws._getattr_unchecked_(attr).initialize_if_not()
        ws._setattr_unchecked_(attr,
                               type(ws._getattr_unchecked_(attr).value)(val))
    elif type(val) in _group_types:
        ws._getattr_unchecked_(attr).initialize_if_not()
        ws._setattr_unchecked_(attr, val)
    else:
        raise RuntimeError(f"'{attr}' is not on Workspace, and "
                           f"{type(val)} is not a Workspace Group")


def _CallbackOperator(func, ins, outs, extras=[], fnstr="Undefined"):
    extras.extend(ins)

    def fn(ws):
        try:
            out = func(*[ws._getattr_unchecked_(x).value for x in ins])

            if len(outs) == 1:
                safe_set(ws, outs[0], out)
            elif len(outs) > 1:
                for i in range(len(outs)):
                    attr = outs[i]
                    val = out[i]
                    safe_set(ws, attr, val)
        except Exception as e:
            raise RuntimeError(
                f"Raised internal error:\n{e}\n\n"
                "The original function (if known) reads:\n"
                f"{fnstr}"
                f"\n\nThe callback is:\nInput: {extras}\nOutput: {outs}"
            )

    return cxx.CallbackOperator(fn, extras, outs)


def _callback_operator_function(funcdef):
    assert isinstance(funcdef, FunctionDef), "Must be a function"
    return [x.arg for x in funcdef.args.args]


def _callback_operator_return(expr):
    value = expr.value

    if isinstance(value, Name):
        return [value.id]
    elif isinstance(value, Tuple) or isinstance(value, List):
        return [x.id if isinstance(x, Name) else _errvar
                for x in value.elts]
    return []


def _find_return(body):
    for expr in body:
        if isinstance(expr, Return):
            return _callback_operator_return(expr)
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

    return None


def _rename_var(name_list, **kwargs):
    """ Renames things in name_list by the kwargs
    """
    if name_list is None:
        return []

    for i in range(len(name_list)):
        n = name_list[i]
        n = kwargs.get(n, n)
        name_list[i] = n
    return name_list


def _callback_operator(func, inputs, outputs, extras, **kwargs):
    """ Internal source code parser
    """
    srccod = getsource(func)
    srccod = unindent(srccod)
    srcast = parse(srccod)

    assert len(srcast.body) == 1

    code_body = srcast.body[0]

    args = (
        inputs
        if len(inputs)
        else _rename_var(_callback_operator_function(code_body), **kwargs)
    )

    retvals = (
        outputs
        if len(outputs)
        else _rename_var(_find_return(code_body.body), **kwargs)
    )

    if _errvar in retvals:
        raise RuntimeError(f"Missing return name in:\n{srccod}\n\n"
                           f"The return-values are {retvals}")

    return _CallbackOperator(func, args, retvals, extras, srccod)


def callback_operator(
    func=None, *, inputs=[], outputs=[], extras=[], **kwargs
):
    """
    Creates a callback operator
    """
    def parser(fn):
        return _callback_operator(fn, inputs, outputs, extras, **kwargs)

    if func is None:
        return parser
    else:
        return _callback_operator(func, inputs, outputs, extras, **kwargs)
