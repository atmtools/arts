import pyarts3.arts as cxx


def recursion(orig, name="pyarts.arts"):
    recr_ = type(cxx.path), type(cxx.Index)
    func = type(cxx.globals.time_report), type(cxx.Index.__init__)

    lst = []

    for m in dir(orig):
        try:
            obj = eval("orig." + m)

            if isinstance(obj, recr_):
                lst.extend(recursion(obj, name + f".{m}"))
            elif isinstance(obj, func):
                for signature in obj.__nb_signature__:
                    if "::" in signature[0]:
                        lst.append(name + f".{m}: {signature[0]}")
        except Exception:
            lst.append(name + f".{m}")
    return lst


errors = recursion(cxx)

for error in errors:
    print(error, end="\n\n")

assert len(errors) == 0, (
    f"Found {len(errors)} signatures with remaining '::' in them.\n\n"
    "This means that there are missing type information when creating the python signature in the C++ bindings.\n\n"
    "The types are not defined in the python interface.\n\n"
    "Please investigate the above and correct these signatures."
)
