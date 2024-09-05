from ast import parse, FunctionDef, Return, Name, Tuple, List, ClassDef
from inspect import getsource
import pyarts.arts as cxx
from pyarts.workspace.utility import unindent as unindent


_errvar = "Unknown Value"


_group_types = [eval(f"cxx.{x}") for x in
                list(cxx.globals.workspace_groups().keys())]


def _CallbackOperator(func, ins, outs, fnstr="Undefined"):
    def fn(ws):
        try:
            out = func(*[ws.get(x) for x in ins])

            if len(outs) == 1:
                setattr(ws, outs[0], out)
            elif len(outs) > 1:
                for i in range(len(outs)):
                    attr = outs[i]
                    val = out[i]
                    setattr(ws, attr, val)
        except Exception as e:
            raise RuntimeError(
                f"Raised internal error:\n{e}\n\n"
                "The original function (if known) reads:\n"
                f"{fnstr}"
                f"\n\nThe callback is:\nInput: {ins}\nOutput: {outs}"
            )
    return cxx.CallbackOperator(fn, ins, outs)


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
    return None


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


def _callback_operator(func, inputs, outputs):
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
        else _callback_operator_function(code_body)
    )

    retvals = (
        outputs
        if len(outputs)
        else _find_return(code_body.body)
    )

    if _errvar in retvals:
        raise RuntimeError(f"Error in:\n{srccod}\n\n"
                           "Missing named return values.\n"
                           "Consider setting `outputs` in decorator.")

    return _CallbackOperator(func, args, retvals, srccod)


def callback_operator(
    func=None, *, inputs=[], outputs=[]
):
    """
    Creates a callback operator
    """
    def parser(fn):
        return _callback_operator(fn, inputs, outputs)

    if func is None:
        return parser
    else:
        return _callback_operator(func, inputs, outputs)
