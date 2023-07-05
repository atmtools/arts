import pyarts
import sys


def sort_ignore_case(entities):
    sorted = [str(e.name) for e in entities]
    sorted.sort(key=lambda x: x.lower())
    return sorted


def agendas():
    vars = sort_ignore_case(pyarts.arts.globals.get_agenda_data())
    txt = """.. hlist::
"""
    for var in vars:
        txt += f"    * :attr:`~pyarts.workspace.Workspace.{var}`\n"
    return txt


def variables():
    vars = sort_ignore_case(pyarts.arts.globals.get_wsv_data())
    txt = """.. hlist::
"""
    for var in vars:
        txt += f"    * :attr:`~pyarts.workspace.Workspace.{var}`\n"
    return txt


def groups():
    vars = sort_ignore_case(pyarts.arts.globals.get_wsv_groups())
    txt = """.. hlist::
"""
    for var in vars:
        txt += f"    * :attr:`~pyarts.arts.{var}`\n"
    return txt


def methods():
    vars = sort_ignore_case(pyarts.arts.globals.get_md_data_raw())
    txt = """.. hlist::
"""
    for var in vars:
        txt += f"    * :func:`~pyarts.workspace.Workspace.{var}`\n"
    return txt


if sys.argv[1] == "Agendas":
    print(agendas())
elif sys.argv[1] == "Variables":
    print(variables())
elif sys.argv[1] == "Groups":
    print(groups())
elif sys.argv[1] == "Methods":
    print(methods())
