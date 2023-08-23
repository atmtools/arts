# Import the module
import pyarts


# Create a workspace
ws = pyarts.Workspace()


"""
All workspace methods can take pure python instance or workspace variables
"""
example_index = pyarts.arts.Index(1)
ws.example_index = example_index
ws.Print(example_index, 0)
ws.Print(ws.example_index, 0)
print(type(ws.example_index), type(example_index))

"""
You can call workspace methods with named or positional arguments or any
combination thereof that python accepts
"""
# 1)  Name some arguments
ws.iy_space_agendaSet(option="CosmicBackground")
# 2)  Name some arguments, use positional for the others
ws.iy_space_agendaSet(ws.iy_space_agenda, option="CosmicBackground")
# 3)  Name all arguments
ws.iy_space_agendaSet(
    iy_space_agenda=ws.iy_space_agenda, option="CosmicBackground"
)
# 4)  Use all positional arguments
ws.iy_space_agendaSet(ws.iy_space_agenda, "CosmicBackground")


"""
All arguments that are marked as input must be provided.  An alternative
way to provide input arguments is to set the corresponding workspace
variable before calling the method.  The following are equivalent:
"""
ws.met_mm_backend = [[1e9, 1e8, 1e7, 1e6]]
ws.f_gridMetMM(met_mm_backend=[[1e9, 1e8, 1e7, 1e6]])
print("A  vector of channel frequencies:\n", ws.f_grid)
ws.met_mm_backend = [[1e9, 1e8, 1e7, 1e6]]
ws.f_gridMetMM()
print("A vector of channel frequencies:\n", ws.f_grid)

"""
Some methods take input arguments that do not have a corresponding
workspace variable.  These arguments may have a default value, or they
may be required.  The example function has some input that are default
and this must be modified when calling the function by setting the argument
manually (by position or, as in the example, by name)
"""
ws.f_gridMetMM(freq_spacing=[5e5])
print("A denser vector of channel frequencies:\n", ws.f_grid)

# TESTING
# AS THIS FILE IS RUN TO TEST ARTS, WE NEED TO CHECK THAT
# THE CONTENT OF THE VARIABLES ARE GOOD.
assert (ws.f_grid == [ 889750000,  890250000,  909750000,  910250000,
                      1089750000, 1090250000, 1109750000, 1110250000]).all()
# END TESTING
