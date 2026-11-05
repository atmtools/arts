import re
from contextlib import contextmanager
from io import StringIO
from pathlib import Path
from types import FunctionType
import pyarts3
import pyarts3.arts as cxx
import sys

global_errors = []


@contextmanager
def rst_output(path):
    """Preserve timestamps so Sphinx can reuse unchanged generated pages."""
    buffer = StringIO()
    yield buffer
    content = buffer.getvalue()
    path = Path(path)
    if not path.exists() or path.read_text(encoding="utf-8") != content:
        path.write_text(content, encoding="utf-8")


module_t = type(cxx)
attr_t = property
class_t = type(cxx.Index)
func_t = type(cxx.Index.fromxml)
method_t = type(cxx.Index.savexml)
enum_t = type(cxx.PType)


def doc(x):
    try:
        return x.__doc__ if x.__doc__ else ""
    except Exception as e:
        raise Exception(f"Error in doc for {x}:\n{e}")


def indent(text, num_spaces):
    res = "\n".join(" " * num_spaces + line for line in text.split("\n"))
    return res.replace(" " * num_spaces + '\n', '\n').rstrip() + "\n"


def short_doc(v, var=None, name=None):
    try:
        x = doc(v).split("\n")
        s = x[0]

        i = 1
        while ("->" in s or len(s) == 0) and i < len(x):
            s = x[i]
            i += 1

        if "->" in s or len(s) == 0:
            attr = v.__name__ if hasattr(v, "__name__") else name
            if not is_operator(attr):
                global_errors.append(
                    f"ERROR [ARTS DOC] No short doc found: {var if var is not None else ""}{'.' if var is not None and attr is not None else ""}{attr if attr is not None else ''}"
                )

        # A signature fallback is code, not prose (e.g. its keyword-only *).
        if "->" in s:
            return f"``{s}``"
        return s
    except Exception as e:
        raise Exception(f"Error in short_doc for v={v}, var={var}, name={name}:\n{e}")


def typesof(v):
    x = doc(v).split(".. :class:")

    classes = None
    for c in x:
        if len(c) < 2 or c[0] != '`':
            continue
        if classes is None:
            classes = c.split("`")[1]
        else:
            classes += " | " + c.split("`")[1]

    if classes is not None:
        classes = re.sub(
            r'(?<![\w.~])([A-Za-z_]\w*)(?![\w.])',
            lambda m: f"~pyarts3.arts.{m[0]}"
            if isinstance(getattr(cxx, m[0], None), type) else m[0],
            classes,
        )
    return classes


# Roots that Sphinx resolves through intersphinx rather than through our own
# objects, so they are trusted without importing them here.
external_roots = {
    "builtins",
    "collections",
    "datetime",
    "numpy",
    "scipy",
    "typing",
}

# Bare names that are valid in a type expression without being ARTS classes.
builtin_names = {
    "None",
    "bool",
    "bytearray",
    "bytes",
    "complex",
    "dict",
    "float",
    "int",
    "list",
    "object",
    "set",
    "str",
    "tuple",
}


def unresolved_type_names(types):
    """Names of a type expression that Sphinx cannot turn into a link.

    typesof() has already spelled every ARTS class out as ``pyarts3.arts.<name>``,
    so a bare name that is left over is either a builtin or a mistake -- and any
    mistake here is a Sphinx warning about an unresolvable reference later on.
    """
    bad = []
    for token in re.findall(r"[A-Za-z_][\w.]*", types):
        if "." in token:
            if token.split(".")[0] in external_roots:
                continue
            obj = pyarts3
            for part in token.split(".")[1:]:
                obj = getattr(obj, part, None)
                if obj is None:
                    break
            if isinstance(obj, type):
                continue
        elif token in builtin_names:
            continue
        elif isinstance(getattr(cxx, token, None), type):
            continue
        bad.append(token)
    return bad


def check_types(owner, name, types):
    """Collect the type-info problems of one attribute as generator errors."""
    if types is None:
        global_errors.append(
            f"{owner}.{name} lacks type-info -- end the docstring with a line of "
            f"the form '.. :class:`<the type, spelled out>`'"
        )
        return "object"
    for bad in unresolved_type_names(types):
        global_errors.append(
            f"{owner}.{name} has type-info that Sphinx cannot resolve: "
            f"{bad} (in {types})"
        )
    return types


def class_ref(t):
    """Reference target for a value's type.

    A generated page sets its currentmodule to the module it documents, which
    for a submodule page is where Sphinx then looks up a bare class name.  The
    classes are documented under pyarts3.arts, so they need their module
    spelled out.  Builtins keep their bare name, which is what the Python
    inventory holds them under.
    """
    mod = getattr(t, "__module__", "") or ""
    return f"{mod}.{t.__name__}" if mod.startswith("pyarts3") else t.__name__


def retypeof(v):
    x = v.split(" | ")

    def fixcur(s):
        return f":class:`{s}`" if s != "None" else ":obj:`None`"

    out = ""
    cur = ""
    first = True
    for c in x:
        if first:
            first = False
        else:
            out += " or "

        cur = ""
        for b in re.split('(\\W+)', c):
            if b == "":
                continue

            if b[0] == "[" or b[0] == ']':
                out += f"{fixcur(cur)}\\{b[0]}"
                cur = b[1:].strip() if len(b) > 1 else ""
                continue

            if "," == b[0]:
                out += f"{fixcur(cur)}, "
                cur = b[1:].strip() if len(b) > 1 else ""
                continue
            cur += b
        if len(cur) > 0:
            out += fixcur(cur)

    return out


def func(name, var):
    try:
        return {"name": name, "short": short_doc(getattr(var, name), var.__name__, name), "types": typesof(getattr(var, name)), "doc": doc(getattr(var, name))}
    except Exception as e:
        raise Exception(f"Error in func for name={name}, var={var}:\n{e}")


def is_internal(name):
    try:
        return not is_operator(name) and (
            name.startswith("_") or name.endswith("__") or ":" in name
        )
    except Exception as e:
        raise Exception(f"Error in is_internal for {name}:\n{e}")


def is_operator(name):
    operators = [
        "__call__",
        "__abs__",
        "__add__",
        "__and__",
        "__contains__",
        "__divmod__",
        "__eq__",
        "__floordiv__",
        "__ge__",
        "__getitem__",
        "__gt__",
        "__iadd__",
        "__iand__",
        "__ifloordiv__",
        "__imatmul__",
        "__imod__",
        "__imul__",
        "__ior__",
        "__ipow__",
        "__isub__",
        "__iter__",
        "__itruediv__",
        "__le__",
        "__len__",
        "__lt__",
        "__matmul__",
        "__mod__",
        "__mul__",
        "__ne__",
        "__hash__",
        "__or__",
        "__pow__",
        "__radd__",
        "__rand__",
        "__rdivmod__",
        "__rfloordiv__",
        "__rmatmul__",
        "__rmod__",
        "__rmul__",
        "__ror__",
        "__rpow__",
        "__rsub__",
        "__rtruediv__",
        "__setitem__",
        "__delitem__",
        "__truediv__",
        "__array__",
        "__init__",
        "__format__",
        "__repr__",
        "__str__",
    ]
    return name in operators


def workspace_page(name):
    # ReadXML and readxml are distinct APIs on case-insensitive filesystems too.
    return "readxml_method" if name == "readxml" else name


def loop_over_class(cls, mod, pure_overview=False):
    try:
        attributes = {}
        functions = {}
        methods = {}
        values = {}
        operators = {}

        for name in dir(cls):
            var = getattr(cls, name)
            if is_internal(name) or (
                isinstance(var, type) and name in {"ItemIterator", "KeyIterator", "ValueIterator"}
            ):
                pass
            elif is_operator(name):
                operators[name] = func(name, cls)
            elif isinstance(var, attr_t):
                attributes[name] = func(name, cls)
            elif isinstance(var, (func_t, FunctionType)):
                functions[name] = func(name, cls)
            elif isinstance(var, method_t):
                methods[name] = func(name, cls)
            else:
                # Nested enums are types too; they get an autoclass block of
                # their own on the parent page, further down.
                values[name] = var

        short = short_doc(cls)

        str = f"{cls.__name__}\n{'#'*len(cls.__name__)}\n\n"

        str += ".. currentmodule:: " + mod + "\n\n"
        str += ".. autoclass:: " + cls.__name__ + "\n\n"

        if (
            len(methods)
            or len(functions)
            or len(attributes)
            or len(values)
            or len(operators)
        ):

            str += f"  .. rubric:: Overview\n\n"
            str += f"  .. list-table::\n\n"

            for n in methods:
                str += f"    * - Method\n"
                str += f"      - :func:`~{mod}.{cls.__name__}.{n}`\n"
                str += f"      - {methods[n]['short']}\n"

            for n in functions:
                str += f"    * - Static Method\n"
                str += f"      - :func:`~{mod}.{cls.__name__}.{n}`\n"
                str += f"      - {functions[n]['short']}\n"

            for n in attributes:
                attributes[n]["types"] = check_types(
                    cls.__name__, n, attributes[n]["types"])
                str += f"    * - {retypeof(attributes[n]["types"])}\n"
                str += f"      - :attr:`~{mod}.{cls.__name__}.{n}`\n"
                str += f"      - {attributes[n]['short']}\n"

            for n in values:
                str += f"    * - Static Data\n"
                str += f"      - ``{mod}.{cls.__name__}.{n}``\n"
                if isinstance(values[n], type):
                    str += f"      - :class:`~{mod}.{cls.__name__}.{n}`\n"
                else:
                    str += f"      - {repr(values[n])} (:class:`~{class_ref(type(values[n]))}`)\n"

            for n in operators:
                str += f"    * - Operator\n"
                str += f"      - :func:`~{mod}.{cls.__name__}.{n}`\n"
                str += f"      - {operators[n]['short']}\n"

        str += "\n"
        str += "  .. rubric:: Constructors\n\n"
        str += "  .. automethod:: __init__\n"
        str += "     :noindex:\n\n"

        if not pure_overview:
            if len(methods):
                str += "  .. rubric:: Methods\n\n"
                for n in methods:
                    str += (
                        f"  .. automethod:: {cls.__name__}.{methods[n]['name']}\n"
                    )
                str += "\n"

            if len(functions):
                str += "  .. rubric:: Static Methods\n\n"
                for n in functions:
                    str += f"  .. automethod:: {cls.__name__}.{functions[n]['name']}\n"
                str += "\n"

            if len(attributes):
                str += "  .. rubric:: Attributes\n\n"
                for n in attributes:
                    str += f"  .. attribute:: {cls.__name__}.{attributes[n]['name']}\n"
                    str += f"     :type: {attributes[n]['types']}\n"
                    str += f"\n"
                    str += f"{indent(attributes[n]['doc'], 5)}\n\n"
                str += "\n"

            if len(operators):
                str += "  .. rubric:: Operators\n\n"
                for n in operators:
                    str += f"  .. automethod:: {cls.__name__}.{operators[n]['name']}\n"
                str += "\n"
        else:
            str += ".. toctree::\n"
            str += "  :hidden:\n\n"

            for n in methods:
                str += f"  {mod}.{cls.__name__}.{workspace_page(n)}\n"

            for n in functions:
                str += f"  {mod}.{cls.__name__}.{workspace_page(n)}\n"

            for n in attributes:
                str += f"  {mod}.{cls.__name__}.{workspace_page(n)}\n"

            for n in operators:
                str += f"  {mod}.{cls.__name__}.{workspace_page(n)}\n"

        for n, value in values.items():
            if isinstance(value, type):
                str += f"\n.. autoclass:: {mod}.{cls.__name__}.{n}\n   :members:\n\n"

        return str, short
    except Exception as e:
        raise Exception(
            f"Error in loop_over_class cls={cls}, mod={mod}, pure_overview={pure_overview}:\n{e}")


def loop_over_module(mod):
    try:
        modules = {}
        classes = {}
        attributes = {}
        functions = {}
        methods = {}
        special = []
        values = {}

        for name in dir(mod):
            var = getattr(mod, name)
            fullname = f"{mod.__name__}.{name}"

            if is_internal(name):
                special.append(fullname)
            elif isinstance(var, module_t):
                modules[fullname] = loop_over_module(var)
            elif isinstance(var, enum_t):
                classes[fullname] = (
                    f"{name}\n{'#' * len(name)}\n\n.. autoclass:: {fullname}\n   :members:\n",
                    (var.__doc__ or "Enumeration.").split("\n")[0],
                )
            elif isinstance(var, class_t):
                classes[fullname] = loop_over_class(var, mod.__name__)
            elif isinstance(var, attr_t):
                attributes[fullname] = func(name, mod)
            elif isinstance(var, (func_t, FunctionType)):
                functions[fullname] = func(name, mod)
            elif isinstance(var, method_t):
                methods[fullname] = func(name, mod)
            elif not isinstance(var, enum_t):
                values[fullname] = var

        return {
            "modules": modules,
            "classes": classes,
            "attributes": attributes,
            "functions": functions,
            "methods": methods,
            "special": special,
            "values": values,
            "short": doc(mod).split("\n")[0] if mod.__doc__ else "",
        }
    except Exception as e:
        raise Exception(f"Error in loop_over_module for mod={mod}:\n{e}")


data = loop_over_module(cxx)


def all_page_names(data, out=None):
    """Every RST page that will be written, to spot name clashes."""
    if out is None:
        out = set()
    for name, sub in data["modules"].items():
        out.add(name)
        all_page_names(sub, out)
    out.update(data["classes"])
    out.update(data["functions"])
    out.update(data["methods"])
    return out


pages = all_page_names(data)


def module_page(name):
    """The page of a module, kept apart from pages that differ only by case.

    A class and a submodule that differ only in case would fight over a single
    file on case-insensitive filesystems, and whichever is written second wins:
    the module page, with the automodule directive that is the only source of
    its py:mod target, is simply gone.  Sphinx takes the target from the
    directive, not from the file name, so the module page can be renamed while
    every :py:mod: reference to it keeps working.
    """
    if any(n != name and n.lower() == name.lower() for n in pages):
        return f"{name}-module"
    return name


def create_func_rst(name, path, mod):
    try:
        with rst_output(f"{path}/{name}.rst") as f:
            f.write(f"{name}\n{'='*len(name)}\n\n")
            f.write(f".. currentmodule:: {mod}\n\n")
            f.write(f".. automethod:: {name}\n\n")
    except Exception as e:
        raise Exception(
            f"Error in create_func_rst for name={name}, path={path}, mod={mod}:\n{e}")


def create_class_rst(data, path, mod):
    try:
        with rst_output(f"{path}/{mod}.rst") as f:
            f.write(data)
    except Exception as e:
        raise Exception(
            f"Error in create_class_rst for data={data}, path={path}, mod={mod}:\n{e}")


def create_rst(data, path, mod, page=None):
    try:
        if page is not None and page != mod:
            # Drop the page this module claimed before its clash with a same-
            # cased class was noticed.  On a case-insensitive filesystem that
            # file is also where the class page would go, and it keeps its old
            # name however often it is rewritten.
            Path(f"{path}/{mod}.rst").unlink(missing_ok=True)

        with rst_output(f"{path}/{page or mod}.rst") as f:
            f.write(f"{mod}\n{'='*len(mod)}\n\n")

            f.write(".. toctree::\n")
            f.write("  :hidden:\n\n")

            for name in data["modules"]:
                f.write(f"  {module_page(name)}\n")
                create_rst(data["modules"][name], path, name, module_page(name))

            for name in data["classes"]:
                f.write(f"  {name}\n")
                create_class_rst(data["classes"][name][0], path, name)

            for name in data["functions"]:
                f.write(f"  {name}\n")
                create_func_rst(name, path, mod)

            for name in data["methods"]:
                f.write(f"  {name}\n")
                create_func_rst(name, path, mod)

            for name in data["attributes"]:
                f.write(f"  {name}\n")

            f.write(f"\n.. automodule:: {mod}\n\n")

            if len(data["modules"]):
                f.write("\n.. rubric:: Modules\n\n")
                f.write("\n.. list-table::\n\n")
                for name in data["modules"]:
                    f.write(f"  * - :py:mod:`~{name}`\n")
                    f.write(f"    - {data['modules'][name]['short']}\n")
                    create_rst(data["modules"][name], path, name, module_page(name))

            if len(data["classes"]):
                f.write("\n.. rubric:: Classes\n\n")
                f.write("\n.. list-table::\n\n")
                for name in data["classes"]:
                    f.write(f"  * - :class:`~{name}`\n")
                    f.write(f"    - {data['classes'][name][1]}\n")
                    create_class_rst(data["classes"][name][0], path, name)

            if len(data["functions"]):
                f.write("\n.. rubric:: Functions\n\n")
                f.write("\n.. list-table:: \n\n")
                for name in data["functions"]:
                    f.write(f"  * - :func:`~{name}`\n")
                    f.write(f"    - {data['functions'][name]["short"]}\n")

            if len(data["methods"]):
                f.write("\n.. rubric:: Methods\n\n")
                f.write("\n.. list-table:: \n\n")
                for name in data["methods"]:
                    f.write(f"  * - :func:`~{name}`\n")
                    f.write(f"    - {data['methods'][name]["short"]}\n")

            if len(data["attributes"]):
                f.write("\n.. rubric:: Attributes\n\n")
                f.write("\n.. list-table:: \n\n")
                for name in data["attributes"]:
                    f.write(f"  * - :attr:`~{name}`\n")
                    f.write(f"    - {data['attributes'][name][1]}\n")

            if len(data["values"]):
                f.write("\n.. rubric:: Static Data\n\n")
                f.write("\n.. list-table::\n\n")
                for n in data["values"]:
                    f.write(f"    * - ``{n}``\n")
                    f.write(
                        f"      - :class:`~{class_ref(type(data['values'][n]))}`\n"
                    )
                    f.write(f"      - {repr(data["values"][n])}\n")
    except Exception as e:
        raise Exception(
            f"Error in create_rst for data={data}, path={path}, mod={mod}:\n{e}")


def create_workspace_rst(path):
    try:
        ws = pyarts3.Workspace

        data = loop_over_class(ws, "pyarts3.workspace", True)
        create_class_rst(data[0], path, f"pyarts3.workspace.Workspace")

        for name in dir(ws):
            if is_internal(name):
                # Remove obsolete pages from earlier generator versions.
                Path(f"{path}/pyarts3.workspace.Workspace.{name}.rst").unlink(missing_ok=True)
                continue
            attr = getattr(ws, name)

            if isinstance(attr, attr_t):
                types = check_types("Workspace", name, typesof(attr))
                with rst_output(
                    f"{path}/pyarts3.workspace.Workspace.{workspace_page(name)}.rst"
                ) as f:
                    f.write(f"{name}\n{'='*len(name)}\n\n")
                    f.write(f".. currentmodule:: pyarts3.workspace\n\n")
                    f.write(f".. attribute:: Workspace.{name}\n")
                    f.write(f"   :type: {types}\n\n")
                    f.write(f"{indent(doc(attr), 3)}\n")
            elif isinstance(attr, (method_t, FunctionType)):
                with rst_output(
                    f"{path}/pyarts3.workspace.Workspace.{workspace_page(name)}.rst"
                ) as f:
                    f.write(f"{name}\n{'='*len(name)}\n\n")
                    f.write(".. currentmodule:: pyarts3.workspace\n\n")
                    f.write(".. automethod:: Workspace." + name + "\n\n")
            elif is_operator(name):
                with rst_output(
                    f"{path}/pyarts3.workspace.Workspace.{workspace_page(name)}.rst"
                ) as f:
                    f.write(f"{name}\n{'='*len(name)}\n\n")
                    f.write(".. currentmodule:: pyarts3.workspace\n\n")
                    f.write(".. automethod:: Workspace." + name + "\n\n")
    except Exception as e:
        raise Exception(f"Error in create_workspace_rst for path={path}:\n{e}")


create_rst(data, sys.argv[1], "pyarts3.arts")

create_workspace_rst(sys.argv[1])

if len(global_errors) != 0:
    print("\n".join(dict.fromkeys(global_errors)))
    print(f"ERRORS FOUND: {len(set(global_errors))}")
    sys.exit(1)
