import os
import sys
from ast import (
    parse,
    Expr,
    Call,
    Name,
    Assign,
    Return,
    Attribute,
    literal_eval,
    FunctionDef,
    Expression,
    unparse,
)
from inspect import getsource, getsourcelines, getfile
from pyarts.workspace.utility import unindent as unindent
from pyarts.arts.globals import workspace_methods, workspace_variables
from pyarts.arts import Agenda, Method
import pyarts.arts as cxx


_group_types = [eval(f"cxx.{x}") for x in list(cxx.globals.workspace_groups())]
_wsvs = cxx.globals.workspace_variables()


class Workspace(cxx._Workspace):
    """
    A wrapper for the C++ workspace object
    """

    def __getattribute__(self, attr):
        if attr.startswith("__"):
            object.__getattribute__(self, attr)

        return super().__getattribute__(attr)

    def __getattr__(self, attr):
        if super().has(attr):
            return super().get(attr)

        raise AttributeError(
            f"'Workspace' object has no attribute '{attr}'")

    def __setattr__(self, attr, value):
        if self.has(attr):
            t = type(self.get(attr))
            if not isinstance(value, t):
                value = t(value)
            self.set(attr, value)
        else:
            if attr in _wsvs:
                super().__setattr__(attr, value)
            elif type(value) in _group_types:
                self.set(attr, value)
            else:
                raise AttributeError(
                    f"'Workspace' object has no attribute '{attr}'")

    def __delattr__(self, attr):
        if attr == '__class__':
            raise AttributeError("You cannot delete __class__")
        if self.has(attr):
            self.set(attr, type(self.get(attr))())
        else:
            raise AttributeError(
                f"'Workspace' object has no attribute '{attr}'")

    def __copy__(self):
        x = super().__copy__()
        out = Workspace()
        out.swap(x)
        return out

    def __deepcopy__(self, *args):
        x = super().__deepcopy__(*args)
        out = Workspace()
        out.swap(x)
        return out

    def __getstate__(self):
        return {"Workspace": super().__getstate__()}

    def __setstate__(self, d):
        super().__setstate__(d["Workspace"])


def _get_attrs(expr):
    """
    Return the name and the value of an attribute
    """
    return expr.value.id, expr.attr


def _eval(expr, state):
    """
    Return a value of an expression
    """
    return eval(
        compile(Expression(body=expr), "<arts_agenda>", "eval"), state
    )


def _assign_parser(expr, ws, state):
    """
    Creates the method that assings a value to the workspace by copy or set

    The current global state is used when RHS is a Name or a Call,
    otherwise literal evaluation is performed to return the value

    Currently set up to handle:

    1) ws.LHS = RHS
    2) ws.LHS = ws.RHS
    3) ws.LHS = <literal>
    4) ws.LHS = <call>

    where LHS is any workspace variable.  Named or not.

    RHS has the following limitations:

    1) RHS must be a value in globals(), its value is parsed and copied.
       The copied value must be possible to assign to LHS.  So if LHS is
       an existing workspace variable, RHS can construct it.  Otherwise,
       RHS must be a workspace group.
    2) None, the assignment will fail upon execution of the method if the two
       workspace variables are of different workspace groups.
    3) <literal> is the catch-all.  If nothing else works, the <literal>
       expression is evaluated at the end to try and evaluate its value before
       returning the method.
    4) Any function call that has a single return value.  The function call
       lives in the globals().  The return value has the same limitation as RHS
       in case 1).
    """
    output_target = expr.targets

    if len(output_target) != 1:
        return "Can only assign to a single target at a time"

    output_target = output_target[0]
    if not isinstance(output_target, Attribute):
        return "Can only assign to attributes"

    output_obj, output_target = _get_attrs(output_target)

    if output_obj != ws:
        return f"Cannot assign to object {output_obj}, did you mean {ws}?"

    value = expr.value

    if isinstance(value, Attribute):
        input_obj, input_target = _get_attrs(value)

        if input_obj == output_obj:
            return Method("Copy", [output_target, input_target])

    tmp = Workspace(False)
    try:
        setattr(tmp, output_target, _eval(value, state))
    except Exception as e:
        return f"Exception: {e}"

    return Method(output_target, getattr(tmp, output_target))


def _method_args(name):
    """
    Return all arguments and their types of a method in order of appearance
    """
    m = workspace_methods().get(name)
    if m is None:
        return f"Unknown method '{m}'"
    out = {}

    for i in m.output:
        out[i] = workspace_variables()[i].type

    for i, t in zip(m.gout, m.gout_type):
        out[i] = t

    for i in m.input:
        if i not in out:
            out[i] = workspace_variables()[i].type

    for i, t in zip(m.gin, m.gin_type):
        out[i] = t

    return out


class _NamedArg:
    """
    Internal type used to name arguments that live on the workspace
    """

    def __init__(self, s):
        self.arg = s


def _call_arg_value(arg, ws, state):
    """
    Get the value of an argument or its name
    """
    if isinstance(arg, Attribute):
        inws, intarget = _get_attrs(arg)
        if inws == ws:
            return _NamedArg(intarget)
    return _eval(arg, state)


def _call_args_parser(call, margs, ws, state):
    """
    Parse all positional arguments of a user-defined method and fill out
    the dict required by the method constructor
    """
    if len(margs) < len(call.args):
        return "Too many inputs"

    out = {}
    for i in range(len(call.args)):
        arg = call.args[i]
        out[margs[i]] = _call_arg_value(arg, ws, state)
    return out


def _call_keywords_parser(kwargs, mdict, call, ws, state):
    """
    Parse and append all named arguments of a user-defined method
    """
    try:
        for keyword in call.keywords:
            arg = keyword.arg

            if arg in kwargs:
                return f'Duplicate argument for "{arg}"'

            if arg not in mdict:
                return f'Unknown argument "{arg}"'

            kwargs[arg] = _call_arg_value(keyword.value, ws, state)

        return kwargs
    except Exception as e:
        return f"Error parsing keyword: {e}"


def _expr_call_parser(call, ws, state):
    """
    Evaluate the call as an ARTS method.  If this is not an ARTS method,
    the CallbackOperator should be copied onto the workspace and the
    CallbackOperatorExecute method should be executed (WIP)
    """
    myws, func = _get_attrs(call.func)

    if myws != ws:
        return f"Bad workspace {myws}, expected {ws}"

    mdict = _method_args(func)

    args = _call_args_parser(call, list(mdict.keys()), ws, state)
    if isinstance(args, str):
        return args

    kwargs = _call_keywords_parser(args, mdict, call, ws, state)
    if isinstance(kwargs, str):
        return kwargs

    methods = []
    call_kwargs = {}
    for k in kwargs:
        kw = kwargs[k]
        t = eval(f"cxx.{mdict[k]}")
        try:
            if isinstance(kw, _NamedArg):
                call_kwargs[k] = f"{kw.arg}"
            else:
                methods.append(Method(f"@{k}", t(kw)))
                call_kwargs[k] = f"@{k}"
        except Exception as e:
            return (
                f"\n{e}\n\n"
                f"Failed to parse {'positional' if k in args else 'named'}"
                f' argument "{k}"'
                f" (arg nr.: {1 + list(kwargs.keys()).index(k)})"
            )

    methods.append(Method(func, [], call_kwargs))
    return methods


def _expr_parser(expr, ws, state):
    """
    An expression is bad unless it is a call, in which case it is forwarded
    """
    if isinstance(expr, Call):
        return _expr_call_parser(expr, ws, state)
    if isinstance(expr, Name):
        return f"Meaningless name.  Did you mean {ws}.{expr.id} = ...?"

    return "Unknown expression value"


def _method_parser(expr, ws, state):
    """
    Parse the agenda method list.  We can only have assignments and expressions
    """
    if isinstance(expr, Expr):
        return _expr_parser(expr.value, ws, state)
    if isinstance(expr, Assign):
        return _assign_parser(expr, ws, state)
    if isinstance(expr, Return):
        return "Cannot return from an Agenda"

    return "Unknown expression encountered parsing method"


def _return_workspace_methods(code_body, ws, state):
    """
    Core workspace method interpreter returning a list of Method and str

    If any str is there, the parsing failed and a future test will produce
    the appropriate error message
    """
    out = []

    for expr in code_body:
        out.append(_method_parser(expr, ws, state))

    return out


def _return_workspace_argname(ast_code):
    """
    Returns the name of the workspace
    """
    args = ast_code.args
    if len(args.args) == 1:
        arg = args.args[0]
        return arg.arg
    return None


def _agenda_or(methods, func, src, fn, lineno):
    """
    Returns an agenda or deal with the error
    """
    agenda = Agenda(func.name)

    for i in range(len(methods)):
        if isinstance(methods[i], Method):
            agenda.add(methods[i])
        elif isinstance(methods[i], list):
            for m in methods[i]:
                agenda.add(m)
        elif isinstance(methods[i], str):
            raise SyntaxError(
                f"In agenda decorator parsing:\n{src}\n\n"
                f'Bad expression "{unparse(func.body[i])}"\n{methods[i]}',
                (fn, lineno + func.body[i].lineno - 1, 0, 0),
            )
        else:
            raise SyntaxError(
                f"In agenda decorator parsing:\n{src}\n\n"
                "Error: Cannot understand code",
                (fn, lineno + func.body[i].lineno - 1, 0, 0),
            )

    return agenda


def _arts_agenda(f, ws, fix):
    """Internal source code parser"""
    srccod = getsource(f)
    srccod = unindent(srccod)
    srcast = parse(srccod)

    assert len(srcast.body) == 1, "Not parsable"

    func = srcast.body[0]
    assert isinstance(func, FunctionDef), "Not a function definition"

    workspace = _return_workspace_argname(func)
    assert workspace is not None, "Must have a workspace variable"

    state = f.__globals__
    methods = _return_workspace_methods(func.body, workspace, state)

    fn = getfile(f)
    ln = getsourcelines(f)[-1]
    agenda = _agenda_or(methods, func, srccod, fn, ln)

    if ws:
        agenda.finalize(fix)
        setattr(ws, agenda.name, agenda)

    return agenda


def arts_agenda(func=None, *, ws=None, fix=False):
    """
    Creates an agenda

    If ws is passed, the agenda is also finalized.  If not, it is still not
    clear that the agenda is valid

    Parameters
    ----------
    func : function
        The function to be turned into an Agenda
    ws : ~pyarts.workspace.Workspace, optional
        The workspace to put this onto after finalization, defaults to None
    fix : bool, optional
        Whether to fix missing input/output in finalization, defaults to False
    """

    def parser(fn):
        return _arts_agenda(fn, ws, fix)

    if func is None:
        return parser
    else:
        return _arts_agenda(func, ws, fix)
