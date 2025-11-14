# Import the module
import pyarts3 as pyarts

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
ws.spectral_rad_space_agendaSet(option="UniformCosmicBackground")

"""
That said, we provide interpreted ways to create agendas manually.
The cosmic background radiation agenda of the method call above
can be created by the following, decorated code:
"""


@pyarts.workspace.arts_agenda
def cosmic_background(ws):
    ws.spectral_radUniformCosmicBackground(
        spectral_rad=ws.spectral_rad,
        freq_grid=ws.freq_grid,
    )
    ws.spectral_rad_jacEmpty()
    ws.Ignore(ws.ray_path_point)


ws.spectral_rad_space_agenda = cosmic_background
print(cosmic_background)
print(ws.spectral_rad_space_agenda)


"""
The local namespace of your python instance will now contain
the variable "cosmic_background".  This is of type Agenda but is unchecked.
The assignment to "ws.iy_space_agenda" will cause the agenda to
be created, set, and checked on the current workspace.  This checking might
fail if, for instance, a required input or output is missing.

It is possible to create a checked agenda that lives on the both the python
and on the pyarts workspace by specifying which workspace it should be set
on as a named argument to the agenda decorator:
"""


@pyarts.workspace.arts_agenda(ws=ws)
def spectral_rad_space_agenda(ws):
    ws.spectral_radUniformCosmicBackground()
    ws.spectral_rad_jacEmpty()
    ws.Ignore(ws.ray_path_point)


print(spectral_rad_space_agenda)
print(ws.spectral_rad_space_agenda)

"""
Lastly, an advanced feature to skip typing ignores is to tell the property to
deal with this for you by fixing the agenda.  Note that this might give you bad
results if the agenda is ever changed in a future update of ARTS as the method
will simply append all new input.  Still, it is convenient:
"""


@pyarts.workspace.arts_agenda(ws=ws, fix=True)
def spectral_rad_space_agenda(ws):
    ws.spectral_radUniformCosmicBackground()
    ws.spectral_rad_jacEmpty()


print(spectral_rad_space_agenda)
