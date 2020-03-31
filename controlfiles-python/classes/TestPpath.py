from pyarts.workspace import Workspace
from pyarts.classes.Ppath import Ppath
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()

# Get the agenda
agenda = from_workspace(ws.ppath)

assert isinstance(agenda, Ppath), "Bad read"
