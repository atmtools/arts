# Import the module
import pyarts
import numpy as np  # For some of the examples


# Create a workspace
ws = pyarts.Workspace()

"""
Indices can be set from the built-in python int.

It is necessary to create an index object first, since the workspace
variable must know the type.
"""
ws.example_index = pyarts.arts.Index(1)
print("Should contain 1:           ", ws.example_index)

"""
Once an index exist on the workspace, you can set it directly from int.
"""
ws.example_index = 2
print("Should contain 2:           ", ws.example_index)

"""
There are several predefined indices in ARTS.  For example, the ybatch_index.
You do not have to create an object for these as the type is already known.
"""
ws.ybatch_index = 3
print("Should contain 3:           ", ws.ybatch_index)


"""
What follows are some examples of how to set other types of workspace
variables.  The general rule is that you can set a workspace variable
from any python object that has a corresponding ARTS type.  For
example, you can set a string from a python string, a matrix from a
numpy array, etc.

All of the variables below are initialization-time known workspace variables,
so there is no need to create the corresponding pyarts class instance first.
"""

# Numeric
ws.g0 = 9.81  # from builtin float
print("Should contain 9.81:        ", ws.g0)

ws.g0 = 10   # from builtin int
print("Should contain 10:          ", ws.g0)

ws.g0 = np.float64(9.81)  # from numpy float
print("Should contain 9.81, again: ", ws.g0)


# String
ws.iy_unit = "dB"  # from builtin str
print("Should contain dB:          ", ws.iy_unit)

ws.iy_unit = b'avc'  # from builtin bytes
print("Should contain avc:         ", ws.iy_unit)


# Vector
ws.f_grid = [1, 2, 3]  # from builtin list
print("Should contain 1 2 3:       ", ws.f_grid)

ws.f_grid = np.array([4, 5, 6])  # from numpy array
print("Should contain 4 5 6:       ", ws.f_grid)


# Matrix
ws.iy = [[1, 2], [3, 4]]  # from builtin list
print("Should contain\n 1 2\n3 4:\n", ws.iy)

ws.iy = np.array([[5, 6], [7, 8]])  # from numpy array
print("Should contain\n 5 6\n7 8:\n", ws.iy)

"""
Many more types are supported.  Please look at the documentation of the
respective Workspace Group for more information on how to initialize them.
"""

# TESTING
# AS THIS FILE IS RUN TO TEST ARTS, WE NEED TO CHECK THAT
# THE CONTENT OF THE VARIABLES ARE GOOD.
assert np.isclose(ws.example_index, 2)
assert np.isclose(ws.ybatch_index, 3)
assert np.isclose(ws.g0, 9.81)
assert ws.iy_unit == "avc"
assert np.allclose(ws.f_grid, [4, 5, 6])
assert np.allclose(ws.iy, [[5, 6], [7, 8]])
# END TESTING
