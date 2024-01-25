# Import the module
import pyarts


# Create a workspace
ws = pyarts.Workspace()


"""
All workspace methods can take pure python instance or workspace variables
"""
example_index = pyarts.arts.Index(1)
ws.example_index = example_index
print(example_index)
print(ws.example_index)
print(type(ws.example_index), type(example_index))

"""
You can call workspace methods with named or positional arguments or any
combination thereof that python accepts
"""
# 1)  Name some arguments
ws.spectral_radiance_background_space_agendaSet(
    option="UniformCosmicBackground"
)
# 2)  Name some arguments, use positional for the others
ws.spectral_radiance_background_space_agendaSet(
    ws.spectral_radiance_background_space_agenda,
    option="UniformCosmicBackground",
)
# 3)  Name all arguments
ws.spectral_radiance_background_space_agendaSet(
    spectral_radiance_background_space_agenda=ws.spectral_radiance_background_space_agenda,
    option="UniformCosmicBackground",
)
# 4)  Use all positional arguments
ws.spectral_radiance_background_space_agendaSet(
    ws.spectral_radiance_background_space_agenda, "UniformCosmicBackground"
)
