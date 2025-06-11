import pyarts.arts as cxx


def recursion(orig, name="pyarts.arts"):
    recr_ = type(cxx.path), type(cxx.Index)
    func = type(cxx.globals.time_report), type(cxx.Index.__init__)

    lst = []

    for m in dir(orig):
        obj = eval("orig." + m)

        if isinstance(obj, recr_):
            lst.extend(recursion(obj, name + f".{m}"))
        elif isinstance(obj, func):
            for signature in obj.__nb_signature__:
                if "std::" in signature[0]:
                    lst.append(name + f".{m}: {signature[0]}")
    return lst


errors = recursion(cxx)

for error in errors:
    print(error, end="\n\n")

assert len(errors) == 0, (
    f"Found {len(errors)} signatures with remaining 'std::' in them.\n\n"
    "This means that there are missing type information when creating the python signature in the C++ bindings.\n\n"
    "There are two possible reasons that we have seen before for this:  1) either the types are not defined "
    "(e.g., not using 'bind_vector', but using the 'std::vector<>' or 'Array<>'), "
    "or 2) there is a missing include (e.g., one of the many '#include <nanobind/stl/XYZ.h>' is missing).  "
    "There could be other reasons we are not aware of for this to happen.\n\n"
    "Please investigate the above and correct these signatures."
)
