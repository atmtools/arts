"""
This example demonstrates how to run a workspace method.
We have shown this in previous examples, but here we will
explain some details.
"""

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
You can call workspace methods with named or unnamed arguments or any combination
thereof that python accepts
"""
ws.iy_space_agendaSet(ws.iy_space_agenda, "CosmicBackground") # Unnamed
ws.iy_space_agendaSet(ws.iy_space_agenda, option="CosmicBackground") # Mixed
ws.iy_space_agendaSet(option="CosmicBackground", iy_space_agenda=ws.iy_space_agenda) # Named


"""
All arguments that the method marks as output can be omitted.  The
following are equivalent:
"""
ws.iy_space_agendaSet(ws.iy_space_agenda, option="CosmicBackground")
ws.iy_space_agendaSet(option="CosmicBackground") # Omitted iy_space_agenda

"""
All arguments that are marked as input must be provided.  An alternative
way to provide input arguments is to set the corresponding workspace
variable before calling the method.  The following are equivalent:
"""
ws.z_surfaceConstantAltitude(lat_grid=[0, 1, 2], lon_grid=[3, 4, 5])
print("A 3-by-3 matrix of zeroes:\n", ws.z_surface)
ws.lat_grid = [0, 1, 2]
ws.lon_grid = [3, 4, 5]
ws.z_surfaceConstantAltitude()
print("A 3-by-3 matrix of zeroes:\n", ws.z_surface)

"""
Some methods take input arguments that do not have a corresponding
workspace variable.  These arguments may have a default value, or they
may be required.  The following is equivalent to the last call of z_surfaceConstantAltitude:
"""
ws.z_surfaceConstantAltitude(altitude=0)
print("A 3-by-3 matrix of zeroes:\n", ws.z_surface)

# TESTING
# AS THIS FILE IS RUN TO TEST ARTS, WE NEED TO CHECK THAT
# THE CONTENT OF THE VARIABLES ARE GOOD.
assert (ws.z_surface.value == 0).all()
# END TESTING
