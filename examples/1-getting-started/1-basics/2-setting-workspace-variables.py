# Import the module
import pyarts3 as pyarts
import numpy as np  # For some of the examples


# Create a workspace
ws = pyarts.Workspace()

"""
Indices can be set from the built-in python int.

It is necessary to create an index object first, since the workspace
variable must know the type.
"""
ws.ind = pyarts.arts.Index(1)
print("Should contain 1:           ", ws.ind)

"""
Once an index exist on the workspace, you can set it directly from int.
"""
ws.ind = 2
print("Should contain 2:           ", ws.ind)


"""
What follows are some examples of how to set other types of workspace
variables.  The general rule is that you can set a workspace variable
from any python object that has a corresponding ARTS type.  For
example, you can set a string from a python string, a matrix from a
numpy array, etc.
"""

# Numeric
ws.num = pyarts.arts.Numeric(9.81)  # from builtin float
print("Should contain 9.81:        ", ws.num)

ws.num = 10   # from builtin int
print("Should contain 10:          ", ws.num)

ws.num = np.float64(9.81)  # from numpy float
print("Should contain 9.81, again: ", ws.num)

# Numeric
ws.str = pyarts.arts.String("hello")  # from builtin str
print("Should contain 'hello':     ", ws.str)

ws.str = b"ARTS"  # from bytes
print("Should contain 'ARTS':      ", ws.str)


# Vector
ws.vec = pyarts.arts.Vector([1, 2, 3])  # from builtin list
print("Should contain 1 2 3:       ", ws.vec)

ws.vec = np.array([4, 5, 6])  # from numpy array
print("Should contain 4 5 6:       ", ws.vec)


# Matrix
ws.mat = pyarts.arts.Matrix([[1, 2], [3, 4]])  # from builtin list
print("Should contain\n 1 2\n3 4:\n", ws.mat)

ws.mat = np.array([[5, 6], [7, 8]])  # from numpy array
print("Should contain\n 5 6\n7 8:\n", ws.mat)

"""
Many more types are supported.  Please look at the documentation of the
respective Workspace Group for more information on how to initialize them.
"""

# TESTING
# AS THIS FILE IS RUN TO TEST ARTS, WE NEED TO CHECK THAT
# THE CONTENT OF THE VARIABLES ARE GOOD.
assert np.isclose(ws.ind, 2)
assert np.isclose(ws.num, 9.81)
assert np.allclose(ws.vec, [4, 5, 6])
assert np.allclose(ws.mat, [[5, 6], [7, 8]])
# END TESTING
