import pyarts.arts as cxx
import sys

module_t = type(cxx)
attr_t = property
class_t = type(cxx.Index)
func_t = type(cxx.Index.fromxml)
method_t = type(cxx.Index.savexml)
enum_t = type(cxx.PType)


def doc(x):
    return x.__doc__ if x.__doc__ else ""

def short_doc(v):
    x = doc(v).split("\n")
    s = x[0]

    i = 1
    while ("->" in s or len(s) == 0) and i < len(x):
        s = x[i]
        i += 1

    if "->" in s or len(s) == 0:
        print("WARNING: [ARTS DOC] No short doc found: ", v.__name__ if hasattr(v, "__name__") else v)

    return s

def func(name, var):
    return {"name": name, "short": short_doc(getattr(var, name))}

def is_internal(name):
    return name.startswith("_") or name.endswith("__") or ":" in name


def loop_over_class(cls, mod):
    attributes = {}
    functions = {}
    methods = {}
    values = {}
    for name in dir(cls):
        var = getattr(cls, name)
        if is_internal(name):
            pass
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

    if len(methods) or len(functions) or len(attributes) or len(values):

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
        str += f"    * - Attribute\n"
        str += f"      - :attr:`~{mod}.{cls.__name__}.{n}`\n"
        str += f"      - {attributes[n]['short']}\n"

      for n in values:
        str += f"    * - Static Data\n"
        str += f"      - ``{mod}.{cls.__name__}.{n}``\n"
        str += f"      - {repr(values[n])} - :class:`~{type(values[n]).__name__}`\n"

    str += "\n"
    str += "  .. rubric:: Constructors\n\n"
    str += "  .. automethod:: __init__\n"
    str += "     :noindex:\n\n"

    if len(methods):
        str += "  .. rubric:: Methods\n\n"
        for n in methods:
            str += f"  .. automethod:: {cls.__name__}.{methods[n]['name']}\n"
        str += "\n"

    if len(functions):
        str += "  .. rubric:: Static Methods\n\n"
        for n in functions:
            str += f"  .. automethod:: {cls.__name__}.{functions[n]['name']}\n"
        str += "\n"

    if len(attributes):
        str += "  .. rubric:: Attributes\n\n"
        for n in attributes:
            str += f"  .. autoattribute:: {cls.__name__}.{attributes[n]['name']}\n"
        str += "\n"

    return str, short


def loop_over_module(mod):
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


data = loop_over_module(cxx)

def create_func_rst(name, path, mod):
    with open(f"{path}/{name}.rst", 'w') as f:
        f.write(f"{name}\n{'='*len(name)}\n\n")
        f.write(f".. currentmodule:: {mod}\n\n")
        f.write(f".. automethod:: {name}\n\n")

def create_class_rst(data, path, mod):
    with open(f"{path}/{mod}.rst", 'w') as f:
        f.write(data)
    

def create_rst(data, path, mod):
    with open(f"{path}/{mod}.rst", 'w') as f:
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
            f.write(f"      - :class:`~{type(data['values'][n]).__name__}`\n")
            f.write(f"      - {repr(data["values"][n])}\n")


create_rst(data, sys.argv[1], "pyarts.arts")
