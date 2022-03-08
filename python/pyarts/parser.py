"""
ARTS controlfile parser
=======================

This module implements a parse for ARTS controlfile. Its implemented
using lark, which greatly simplifies the parsing. Functions are provided
to transform the parsed controlfile to a Python script.
"""
import numpy as np
import re
from textwrap import indent

from lark import Lark, Transformer, Token

import pyarts.workspace.global_data as global_data


workspace_methods = global_data.get_raw_method_map()
workspace_variables = global_data.get_variables_map()
group_names =  list(global_data.cxx.get_wsv_group_names())

grammar = r"""
    controlfile : statement*

    statement : agenda
              | comment
              | function
              | include
              | agenda_definition
              | agenda_append

    include : "INCLUDE " STRING

    agenda : CNAME "{" statement* "}"

    agenda_definition : "AgendaSet" "(" CNAME ")" "{" statement* "}"

    agenda_append : "ArrayOfAgendaAppend" "(" CNAME ")" "{" statement* "}"

    comment : /[ \t]*#.*\n/

    function : CNAME ("(" arguments  ")")? 

    arguments : named_arguments
              | positional_arguments

    named_arguments : [argument_pair ("," argument_pair? )*]

    argument_pair : comment* CNAME comment*  "=" comment* value comment*

    positional_arguments: (comment* value comment* | comment+) (("," comment*  value)  | comment+)*

    list : "[" ((comment* value comment*)? (comment* ("," ) comment* value | comment+)*
               | nested_list (";" nested_list ";"? comment? | comment)+) "]"

    nested_list : (comment* value comment*) (comment* (",") comment* value | comment)*

    matrix : "[" ( | comment)? (comment? ("," | ";") comment? value | comment)* "]"

    empty_list : "[" "]"

    ?value : STRING
           | SIGNED_FLOAT -> number
           | SIGNED_INT
           | list
           | CNAME
           | comment

    DOUBLE_QUOTED_STRING  : /"[^"]*"/
    SINGLE_QUOTED_STRING  : /'[^']*'/
    STRING : (SINGLE_QUOTED_STRING | DOUBLE_QUOTED_STRING)

%import common.SIGNED_FLOAT
%import common.SIGNED_INT
%import common.CNAME
%import common.WS
%import common.WS_INLINE
%import common.NEWLINE
%ignore WS
"""
arts_parser = Lark(grammar, start="controlfile", debug=False, parser="earley")

################################################################################
# Python representation of syntax elements
################################################################################

replace_array = re.compile(r'array\(([^\)]*)\)')
replace_dtype = re.compile(r'dtype=([^\s]+)')

def to_python(obj, workspace):
    """
    Generic function to write elements of a controlfile AST in
    Python syntax. For classes defined below the to_python member
    function is called. Arrays are printed so that they are parsed as
    numpy arrays. Strings are escaped. Other objects (int, list of int)
    are just converted to a string.

    Arguments:
        obj: Element of controlfile AST to write in Python syntax
        workspace: Variable name to use for workspace.
    """
    if hasattr(obj, "to_python"):
        return obj.to_python(workspace)
    elif isinstance(obj, np.ndarray):
        if obj.size == 0:
            return "[]"
        s = repr(obj)
        s = replace_array.sub(r"np.array(\1)", s)
        s = replace_dtype.sub(r"dtype=np.\1", s)
        return  s
    elif isinstance(obj, str):
        return "\"" + str(obj) + "\""
    else:
        return str(obj)


class WSMCall:
    """
    Represents a call of a WSM.

    Attributes:
        name: Name of the WSM that is called
        args: Positional arguments of the WSM call
        kwargs: Named arguments of the call
    """
    def __init__(self, name, args, kwargs):

        if not name in workspace_methods:
            raise Exception("{} is not a known workspace method.".format(name))

        self.wsm = workspace_methods[name]

        self.wsm_outs = [global_data.get_variable_name(m) for m in self.wsm.outs]
        self.wsm_gouts = list(self.wsm.g_out)
        self.wsm_ins = [global_data.get_variable_name(m) for m in self.wsm.ins \
                        if not m in self.wsm.outs]
        self.wsm_gins = list(self.wsm.g_in)
        self.arg_names = self.wsm_outs + self.wsm_gouts + self.wsm_ins + self.wsm_gins

        self.name = name
        self.args = args
        self.kwargs = kwargs

        if not kwargs is None:
            if "in" in kwargs:
                self.kwargs_to_args()

    def kwargs_to_args(self):
        """
        Convert function call from named arguments to positional arguments.
        """
        if self.kwargs == None:
            return None

        args = []

        for n in self.wsm_outs:
            if not n in self.kwargs:
                args.append(WSV(n))
            else:
                args.append(self.kwargs[n])

        for n in self.wsm_gouts:
            args.append(self.kwargs[n])

        for n in self.wsm_ins:
            if not n in self.kwargs:
                args.append(WSV(n))
            else:
                args.append(self.kwargs[n])

        for n in self.wsm_gins:
            if not n in self.kwargs:
                i = self.wsm.g_in.index(n)
                args.append(global_data.convert(self.wsm.g_in_types[i], self.wsm.g_in_default[i]))
            else:
                args.append(self.kwargs[n])


        self.kwargs = None
        self.args = args

    def __repr__(self):
        """
        Print WSM call in ARTS script.
        """
        if self.args is None and self.kwargs is None:
            return self.name + "()\n"

        s = self.name + "("
        if self.kwargs is None:
            for a in self.args[:-1]:
                s += str(a) + ", "
            s += str(self.args[-1]) + ")\n"
        if self.args is None:
            for k in list(self.kwargs.keys())[:-1]:
                s += str(k) + "=" + str(self.kwargs[k]) + ", "
            k = list(self.kwargs.keys())[-1]
            s += str(k) + "=" + str(self.kwargs[k]) + ")\n"
        return s

    def convert_argument(self, name, value):
        """
        Tries to infer type of argument based on types of input and
        generic input.
        """
        if isinstance(value, WSV):
            return value

        if name in self.wsm_ins:
            v = workspace_variables[name]
            value_converted = global_data.convert(v.group, value)
            if not value_converted is None:
                value = value_converted

        if name in self.wsm_gins:
            if len(self.wsm.g_in_types) == 1:
                g = group_names[self.wsm.g_in_types[0]]

                value_converted = global_data.convert(g, value)
                if not value_converted is None:
                    value = value_converted
        return value

    def to_python(self, workspace = "ws"):
        """
        Rewrite function call in Python.
        """
        s = workspace + "." + self.name
        if self.args is None and self.kwargs is None:
            return s + "()\n"
        else:
            s += "("

        if len(self.name) > 6 and self.name[-6:] == "Create":
            if not self.args is None:
                self.args[0] = self.args[0].name
            if not self.kwargs is None:
                k = self.kwargs.keys()
                self.kwargs[k] = self.kwargs[k].name

        if self.kwargs is None:
            for a, n in zip(self.args[:-1], self.arg_names):
                if not isinstance(a, WSV):
                    a = self.convert_argument(n, a)
                s += to_python(a, workspace) + ", "

            a = self.args[-1]
            n = self.arg_names[len(self.args)-1]
            if not isinstance(a, WSV):
                a = self.convert_argument(n, a)
            s += to_python(a, workspace) + ")\n"

        if self.args is None:
            keys = list(self.kwargs.keys())
            for k in keys[:-1]:
                a = self.kwargs[k]
                if not isinstance(a, WSV):
                    a = self.convert_argument(k, a)
                s += str(k) + "=" + to_python(a, workspace) + ", "
            k = keys[-1]
            a = self.kwargs[keys[-1]]
            if not isinstance(a, WSV):
                a = self.convert_argument(k, a)
            s += str(k) + "=" + to_python(a, workspace) + ")\n"
        return s


class AgendaDefinition:
    """
    An agenda defined in a controlfile.

    Attributes:
        name: Name of the agenda
        content: List of statements in the agenda
    """
    def __init__(self, name, content):
        self.name = name
        self.content = content

    def __repr__(self):
        """
        Print agenda definition in controlfile syntax.
        """
        s = "AgendaSet(" + self.name + ") {\n"
        for c in self.content:
            s += str(c)
        s += "}\n"
        return s

    def to_python(self, workspace):
        """
        Print agenda definition in Python syntax.
        """
        s = f"@arts_agenda(ws={workspace})\ndef " + self.name + "({}):\n".format(workspace)
        cs = ""
        for c in self.content:
            cs += to_python(c, workspace)
        s = s + indent(cs, " " * 4)
        s += workspace + "." + self.name + " = " + self.name + "\n\n"
        return s

class AgendaAppend:
    """
    Weird ARTS syntax feature to append agenda to array.

    Attributes:
        name: Name of the agenda to append to
        content: List of statements in the agenda
    """
    def __init__(self, name, content):
        self.name = name
        self.content = content

    def __repr__(self):
        """
        Print agenda definition in controlfile syntax.
        """
        s = "ArrayOfAgendaAppend" + self.name + ") {\n"
        for c in self.content:
            s += str(c)
        s += "}\n"
        return s

    def to_python(self, workspace):
        """
        Print agenda definition in Python syntax.
        """
        s = f"@arts_agenda(ws={workspace})\ndef " + self.name + "({}):\n".format(workspace)
        cs = ""
        for c in self.content:
            cs += to_python(c, workspace)
        s = s + indent(cs, " " * 4) + "\n"
        s += (workspace + ".Append(" + workspace + "." + self.name
              + ", " + self.name + ")\n\n")
        return s

class Comment:
    """
    A comment
    """
    def __init__(self, text):
        self.text = text

    def __repr__(self):
        """
        Print comment in controlfile syntax.
        """
        return str(self.text)

    def to_python(self, workspace):
        """
        Print comment in Python syntax.
        """
        return self.__repr__()


class Include:
    """
    A INCLUDE statement.
    """
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        """
        Print INCLUDE statement in controlfile syntax.
        """
        return "INCLUDE " + "\"" + str(self.name) + "\"\n"

    def to_python(self, workspace):
        """
        Print INCLUDE statement in Python syntax.
        """
        s = "ws.execute_controlfile(\"" + self.name + "\")\n"
        return s

class WSV:
    """
    A workspace variable.

    Attributes:
        name: Name of the WSV
    """
    def __init__(self, name):
        if name in workspace_methods:
            name = camel_to_snake(name)
        self.name = name

    def __repr__(self):
        """
        Print WSV in controlfile syntax.
        """
        return self.name

    def to_python(self, workspace):
        """
        Print WSV in controlfile syntax.
        """
        return workspace + "." + self.name

class Agenda:
    """
    Class to represent the ARTS2 agenda which is the main
    part of a controlfile.

    Attributes:
        name: Name of the agenda
        content: The statement within the agenda.
    """
    def __init__(self, name, content):
        self.name = name
        self.content = content

    def __repr__(self):
        """
        Print agenda in controlfile syntax.
        """
        s = self.name + " {\n"
        for c in self.content:
            s += str(c)
        s += "\n}"
        return s

    def to_python(self, workspace):
        """
        Print agenda in Python syntax.

        """
        s = """
import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda
{} = Workspace(verbosity=0)
""".format(workspace)

        for c in self.content:
            s += to_python(c, workspace)
        return s

class Controlfile:
    """
    Class to represent a whole parse controlfile. Consists
    of a number of comments and one Agenda object.

    Attributes:
        name: Name of the agenda
        content: The statement within the agenda.
    """
    def __init__(self, content):
        self.content = content

    def __repr__(self):
        s = ""
        for c in self.content:
            s += str(c)
        return s

    def to_python(self, workspace):
        s = ""
        for c in self.content:
            s += to_python(c, workspace)
        return s

class ArtsTransformer(Transformer):
    """
    Transformer to transformt lark AST into symbolic representation
    of ARTS controlfile.
    """
    def comment(self, s):
        return Comment(s[0])

    def SIGNED_INT(self, i):
        i = int(i)
        return int(i)

    def include(self, i):
        return Include(i[0])

    def number(self, i):
        return float(i[0])

    def STRING(self, s):
        return s[1:-1]

    def arguments(self, c):
        c = [e for e in c if not isinstance(e, Comment)]
        return c

    def list(self, c):
        c = [e for e in c if not isinstance(e, Comment)]
        return c

    def nested_list(self, c):
        c = [e for e in c if not isinstance(e, Comment)]
        return c

    def function(self, f):
        if len(f) == 1:
            return WSMCall(f[0], None, None)
        else:
            if type(f[1][0]) is list:
                return WSMCall(f[0], f[1][0], None)
            if type(f[1][0]) is dict:
                return WSMCall(f[0], None, f[1][0])

    def positional_arguments(self, c):
        if not type(c) == list:
            c = [c]

        cs = []
        for a in c:
            if isinstance(a, Token) and a.type == "CNAME":
                a = WSV(a)
            cs += [a]

        cs = [e for e in cs if not isinstance(e, Comment)]
        return cs

    def named_arguments(self, c):
        return dict(c)

    def argument_pair(self, p):
        """
        Arguments pairs have a string on the left and a python
        literal or a variable name on the right.
        """
        p = [e for e in p if not isinstance(e, Comment)]
        pl, pr = p
        pl = str(pl)
        if isinstance(pr, Token) and pr.type == "CNAME":
            pr = WSV(pr)
        return (pl, pr)

    def statement(self, s):
        return s[0]

    def controlfile(self, c):
        return Controlfile(c)

    def agenda(self, c):
        return Agenda(c[0], c[1:])

    def agenda_definition(self, a):
        return AgendaDefinition(a[0], a[1:])

    def agenda_append(self, a):
        return AgendaAppend(a[0], a[1:])


################################################################################
# Functions to convert ARTS controlfile to Python
################################################################################

def convert_to_python(controlfile, output, workspace = "ws"):
    with open(controlfile) as f:
        source = f.read()
    tree = arts_parser.parse(source)
    t = ArtsTransformer().transform(tree)
    s = t.to_python(workspace)

    with open(output, "w") as f:
        f.write(s)

pattern = re.compile(r'(?<!^)(?=[A-Z])')
def camel_to_snake(s):
    s = pattern.sub('_', s).lower()
    return s
