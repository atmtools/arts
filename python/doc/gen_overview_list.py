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
    extra = list(pyarts.arts.globals.workspace_agendas_extra().keys())
    main.extend(extra)
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
