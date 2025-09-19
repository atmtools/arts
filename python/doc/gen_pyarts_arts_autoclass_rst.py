import re
import pyarts3
import pyarts3.arts as cxx
import sys

global_errors = []

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

    return classes


def retypeof(v):
    x = v.split(" | ")

    def fixcur(s):
        return f":class:`{s}`" if s != "None" else ":attr:`None`"

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
        "__getstate__",
        "__setstate__",
        "__init__",
        "__format__",
        "__repr__",
        "__str__",
    ]
    return name in operators


def loop_over_class(cls, mod, pure_overview=False):
    try:
        attributes = {}
        functions = {}
        methods = {}
        values = {}
        operators = {}

        for name in dir(cls):
            var = getattr(cls, name)
            if is_internal(name):
                pass
            elif is_operator(name):
                operators[name] = func(name, cls)
            elif isinstance(var, attr_t):
                attributes[name] = func(name, cls)
            elif isinstance(var, func_t):
                functions[name] = func(name, cls)
            elif isinstance(var, method_t):
                methods[name] = func(name, cls)
            elif not isinstance(var, enum_t):
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
                assert attributes[n]["types"] is not None, f"{cls.__name__}.{attributes[n]['name']} lacks type-info - add '.. :class:`class-information`' to a newline of the docstring"
                str += f"    * - {retypeof(attributes[n]["types"])}\n"
                str += f"      - :attr:`~{mod}.{cls.__name__}.{n}`\n"
                str += f"      - {attributes[n]['short']}\n"

            for n in values:
                str += f"    * - Static Data\n"
                str += f"      - ``{mod}.{cls.__name__}.{n}``\n"
                str += f"      - {repr(values[n])} (:class:`~{type(values[n]).__name__}`)\n"

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
                    if attributes[n]['types'] is None:
                        global_errors.append(
                            f"{cls.__name__}.{attributes[n]['name']} lacks type-info - add '.. :class:`class-information`' to a newline of the docstring")
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
                str += f"  {mod}.{cls.__name__}.{n}\n"

            for n in functions:
                str += f"  {mod}.{cls.__name__}.{n}\n"

            for n in attributes:
                str += f"  {mod}.{cls.__name__}.{n}\n"

            for n in operators:
                str += f"  {mod}.{cls.__name__}.{n}\n"

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
            elif isinstance(var, class_t):
                classes[fullname] = loop_over_class(var, mod.__name__)
            elif isinstance(var, attr_t):
                attributes[fullname] = func(name, mod)
            elif isinstance(var, func_t):
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


def create_func_rst(name, path, mod):
    try:
        with open(f"{path}/{name}.rst", "w") as f:
            f.write(f"{name}\n{'='*len(name)}\n\n")
            f.write(f".. currentmodule:: {mod}\n\n")
            f.write(f".. automethod:: {name}\n\n")
    except Exception as e:
        raise Exception(
            f"Error in create_func_rst for name={name}, path={path}, mod={mod}:\n{e}")


def create_class_rst(data, path, mod):
    try:
        with open(f"{path}/{mod}.rst", "w") as f:
            f.write(data)
    except Exception as e:
        raise Exception(
            f"Error in create_class_rst for data={data}, path={path}, mod={mod}:\n{e}")


def create_rst(data, path, mod):
    try:
        with open(f"{path}/{mod}.rst", "w") as f:
            f.write(f"{mod}\n{'='*len(mod)}\n\n")

            f.write(".. toctree::\n")
            f.write("  :hidden:\n\n")

            for name in data["modules"]:
                f.write(f"  {name}\n")
                create_rst(data["modules"][name], path, name)

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
                    create_rst(data["modules"][name], path, name)

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
                        f"      - :class:`~{type(data['values'][n]).__name__}`\n"
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
            attr = getattr(ws, name)

            if isinstance(attr, attr_t):
                with open(
                    f"{path}/pyarts3.workspace.Workspace.{name}.rst", "w"
                ) as f:
                    f.write(f"{name}\n{'='*len(name)}\n\n")
                    f.write(f".. currentmodule:: pyarts3.workspace\n\n")
                    f.write(f".. attribute:: Workspace.{name}\n")
                    f.write(f"   :type: {typesof(attr)}\n\n")
                    f.write(f"{indent(doc(attr), 3)}\n")
            elif isinstance(attr, method_t):
                with open(
                    f"{path}/pyarts3.workspace.Workspace.{name}.rst", "w"
                ) as f:
                    f.write(f"{name}\n{'='*len(name)}\n\n")
                    f.write(".. currentmodule:: pyarts3.workspace\n\n")
                    f.write(".. automethod:: Workspace." + name + "\n\n")
            elif is_operator(name):
                with open(
                    f"{path}/pyarts3.workspace.Workspace.{name}.rst", "w"
                ) as f:
                    f.write(f"{name}\n{'='*len(name)}\n\n")
                    f.write(".. currentmodule:: pyarts3.workspace\n\n")
                    f.write(".. automethod:: Workspace." + name + "\n\n")
    except Exception as e:
        raise Exception(f"Error in create_workspace_rst for path={path}:\n{e}")


create_rst(data, sys.argv[1], "pyarts3.arts")

create_workspace_rst(sys.argv[1])

if len(global_errors) != 0:
    print("\n".join(global_errors))
    print(f"ERRORS FOUND: {len(global_errors)}")
    sys.exit(1)
