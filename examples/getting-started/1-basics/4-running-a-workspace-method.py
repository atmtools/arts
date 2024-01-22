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


"""
All arguments that are marked as input must be provided.  An alternative
way to provide input arguments is to set the corresponding workspace
variable before calling the method.  The following are equivalent:
"""
f = 1e9
# 1
ws.backend_channel_responseGaussian(f_backend=[f], fwhm=[f])
print("Var:", ws.backend_channel_response)
# 2
ws.f_backend = [f]
ws.backend_channel_responseGaussian(fwhm=[f])
print("Var:", ws.backend_channel_response)

"""
Some methods take input arguments that do not have a corresponding
workspace variable.  These arguments may have a default value, or they
may be required.  The example function has one input that is requires and some
that are not. We can change both (by position or, as in the example, by name):
"""
ws.backend_channel_responseGaussian(
    fwhm=[f], grid_width=f / 2.0, grid_npoints=41
)
print("Var:", ws.backend_channel_response)

# TESTING
# AS THIS FILE IS RUN TO TEST ARTS, WE NEED TO CHECK THAT
# THE CONTENT OF THE VARIABLES ARE GOOD.
assert (
    ws.backend_channel_response[0].grids[0]
    == [
        -2.5e08,
        -2.375e08,
        -2.25e08,
        -2.125e08,
        -2e08,
        -1.875e08,
        -1.75e08,
        -1.625e08,
        -1.5e08,
        -1.375e08,
        -1.25e08,
        -1.125e08,
        -1e08,
        -8.75e07,
        -7.5e07,
        -6.25e07,
        -5e07,
        -3.75e07,
        -2.5e07,
        -1.25e07,
        0,
        1.25e07,
        2.5e07,
        3.75e07,
        5e07,
        6.25e07,
        7.5e07,
        8.75e07,
        1e08,
        1.125e08,
        1.25e08,
        1.375e08,
        1.5e08,
        1.625e08,
        1.75e08,
        1.875e08,
        2e08,
        2.125e08,
        2.25e08,
        2.375e08,
        2.5e08,
    ]
).all()
# END TESTING
