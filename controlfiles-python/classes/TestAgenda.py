from pyarts.workspace import Workspace
from pyarts.classes.Agenda import Agenda
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()

# Get the agenda
agenda = from_workspace(ws.abs_xsec_agenda)

assert isinstance(agenda, Agenda), "Bad read"
