import pyarts3 as pyarts
import sys


def sort_ignore_case(entities):
    sorted = [str(e) for e in entities]
    sorted.sort(key=lambda x: x.lower())
    return sorted


def _hlist_num_cols(v):
    return 1 if len(v) < 5 else 2


def agendas():
    existing = pyarts.workspace.Workspace().__dir__()

    main = list(pyarts.arts.globals.workspace_agendas().keys())
    vars = sort_ignore_case(main)
    txt = f""".. hlist::
    :columns: {_hlist_num_cols(vars)}

"""
    for var in vars:
        if var in existing:
            txt += f"    * :attr:`~pyarts3.workspace.Workspace.{var}`\n"
    return txt


def variables():
    existing = pyarts.workspace.Workspace().__dir__()

    vars = sort_ignore_case(list(pyarts.arts.globals.workspace_variables().keys()))
    txt = f""".. hlist::
    :columns: {_hlist_num_cols(vars)}

"""
    for var in vars:
        if var in existing:
            txt += f"    * :attr:`~pyarts3.workspace.Workspace.{var}`\n"
    return txt


def groups():
    existing = pyarts.arts.__dir__()

    vars = sort_ignore_case(list(pyarts.arts.globals.workspace_groups().keys()))
    txt = f""".. hlist::
    :columns: {_hlist_num_cols(vars)}

"""
    for var in vars:
        if var in existing:
            txt += f"    * :class:`~pyarts3.arts.{var}`\n"
    return txt


def _wsv_list(names, existing):
    """Bullet the given workspace variables, or say that there are none."""
    names = [n for n in sort_ignore_case(names) if n in existing]
    if not names:
        return "*None.*\n"

    txt = f""".. hlist::
    :columns: {_hlist_num_cols(names)}

"""
    for name in names:
        txt += f"    * :attr:`~pyarts3.workspace.Workspace.{name}`\n"
    return txt


def dimensions():
    """One section per dimension, listing the variables that share it.

    The label written for each dimension is what variable_dimension_docs() in
    src/workspace_dimensions.cpp links its effective shapes to, so the two have
    to be changed together.
    """
    existing = pyarts.workspace.Workspace().__dir__()
    wsvs = pyarts.arts.globals.workspace_variables()

    txt = ""
    for dim in sort_ignore_case(pyarts.arts.globals.workspace_dimensions().keys()):
        desc = pyarts.arts.globals.workspace_dimensions()[dim].desc

        outer = [n for n, v in wsvs.items() if dim in v.dims]
        inner = [n for n, v in wsvs.items() if dim in v.inner_dims]

        txt += f"""
.. _wsd-{dim}:

{dim}
{'-' * len(dim)}

The {desc}.

.. rubric:: Variables of this size

{_wsv_list(outer, existing)}
.. rubric:: Containers whose elements are of this size

{_wsv_list(inner, existing)}"""
    return txt


def methods():
    existing = pyarts.workspace.Workspace().__dir__()

    vars = sort_ignore_case(list(pyarts.arts.globals.workspace_methods().keys()))
    txt = f""".. hlist::
    :columns: {_hlist_num_cols(vars)}

"""
    for var in vars:
        if var in existing:
            txt += f"    * :func:`~pyarts3.workspace.Workspace.{var}`\n"
    return txt


if sys.argv[1] == "Agendas":
    print(agendas())
elif sys.argv[1] == "Variables":
    print(variables())
elif sys.argv[1] == "Groups":
    print(groups())
elif sys.argv[1] == "Methods":
    print(methods())
elif sys.argv[1] == "Dimensions":
    print(dimensions())
