# Import the module
import pyarts


# Create a workspace
ws = pyarts.Workspace()


"""
Agendas are a collection of functions that define the simulation
of some core component.  For example, the iy_space_agenda defines
how the radiative transfer equation receives radiation from a 
background of space --- commonly just the cosmic background radiation.

These Agendas are defined from a set of common inputs and outputs.

Agendas must consist of methods that combined covers the inputs and
outputs of the agenda itself.  While it is possible to create an
agenda manually, it is much more preferable to use the workspace
methods that exist to create the agenda for you.

Before showing you how to creating your own iy_space_agenda,
this method call can do it for you, and it will probably do
it faster and safer than any manual approach:
"""
ws.iy_space_agendaSet(option="CosmicBackground")

"""
That said, we provide interpreted ways to create agendas manually.
The cosmic background radiation agenda of the method call above
can be created by the following, decorated code:
"""
@pyarts.workspace.arts_agenda
def cosmic_background(ws):
    ws.MatrixCBR(output=ws.iy, f=ws.f_grid)
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)
ws.iy_space_agenda = cosmic_background
print("cosmic_background:", type(cosmic_background))
print("ws.iy_space_agenda:", type(ws.iy_space_agenda.value))

"""
The local namespace of your python instance will now contain
the variable "cosmic_background".  This is of type DelayedAgenda.
The assignment to "ws.iy_space_agenda" will cause the agenda to
be created and set on the workspace.  DelayedAgenda:s are useful
when you want to set up your agendas independently of the workspace.

It is possible to create an actual Agenda object by naming the
workspace input as a named-argument to the decorator:
"""
@pyarts.workspace.arts_agenda(ws=ws)
def cosmic_background_instant(ws):
    ws.MatrixCBR(output=ws.iy, f=ws.f_grid)
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)
ws.iy_space_agenda = cosmic_background_instant
print("cosmic_background_instant:", type(cosmic_background_instant))

"""
The "cosmic_background_instant" variable is now an Agenda instance and
can be used only with the named workspace instance.

You can also just set the agenda to the workspace directly by providing
the "set_agenda=True" argument to the decorator:
"""
@pyarts.workspace.arts_agenda(ws=ws, set_agenda=True)
def cosmic_background_set(ws):
    ws.MatrixCBR(output=ws.iy, f=ws.f_grid)
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)
ws.iy_space_agenda = cosmic_background_set
print("ws.cosmic_background_set:", type(ws.cosmic_background_set.value))

"""
Now the type of the "cosmic_background_set" variable is actually a
WorkspaceVariable, it lives already on the workspace and it has
defined "ws.cosmic_background_set".  Note that this is not very
useful as Agenda:s must be named properly to work in Arts.  Their
names themselves carry a lot of meaning for the internal controlflow.


In fact, it is better to simply name the decorated function
"iy_space_agenda" as this will cause the agenda to be set on
the workspace automatically:
"""
@pyarts.workspace.arts_agenda(ws=ws, set_agenda=True)
def iy_space_agenda(ws):
    ws.MatrixCBR(output=ws.iy, f=ws.f_grid)
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)

"""
Finally, it is not possible by default to create an agenda that
calls pure python functions.  This is because the agenda must
know its inputs and outputs.  If you really want an agenda that
can call pure python functions, you can do it by providing the
"allow_callbacks=True" argument to the decorator:
"""
@pyarts.workspace.arts_agenda(ws=ws, set_agenda=True, allow_callbacks=True)
def iy_space_agenda(ws):
    print("Running iy_space_agenda")
    ws.MatrixCBR(output=ws.iy, f=ws.f_grid)
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)

"""
Note that this is not recommended, as it is not possible to isolate inputs
from outputs when callbacks are involved.  This means that the agenda
above may leak any changes it does to the workspace to the outside world.

The whole point of the Agenda system is to control the inputs and outputs
of the simulation, so it is not recommended to make use of this feature
as it is too powerful in scope.
"""
