import pyarts
import sys


def agendas():
    vars = pyarts.arts.globals.get_agenda_data()
    txt = """.. hlist::
"""
    for var in vars:
        txt += f"    * :attr:`~pyarts.workspace.Workspace.{var.name}`\n"
    return txt


def variables():
    vars = pyarts.arts.globals.get_wsv_data()
    txt = """.. hlist::
"""
    for var in vars:
        txt += f"    * :attr:`~pyarts.workspace.Workspace.{var.name}`\n"
    return txt


def groups():
    vars = pyarts.arts.globals.get_wsv_groups()
    txt = """.. hlist::
"""
    for var in vars:
        txt += f"    * :attr:`~pyarts.arts.{var.name}`\n"
    return txt


def methods():
    vars = pyarts.arts.globals.get_md_data_raw()
    txt = """.. hlist::
"""
    for var in vars:
        txt += f"    * :func:`~pyarts.workspace.Workspace.{var.name}`\n"
    return txt


if sys.argv[1] == "Agendas":
    print(agendas())
elif sys.argv[1] == "Variables":
    print(variables())
elif sys.argv[1] == "Groups":
    print(groups())
elif sys.argv[1] == "Methods":
    print(methods())
