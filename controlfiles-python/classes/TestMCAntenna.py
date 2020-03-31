from pyarts.workspace import Workspace
from pyarts.classes.Vector import Vector
from pyarts.classes.Matrix import Matrix
from pyarts.classes.MCAntenna import MCAntenna
from pyarts.classes import from_workspace

ws = Workspace()

mca = MCAntenna()

mca.type = 1
mca.type = 2
mca.type = 3

mca.sigma_aa = 4.5
mca.sigma_za = 4.5

mca.aa_grid = Vector([1, 2, 3, 4])
mca.za_grid = Vector([1, 2, 3, 4, 5])
mca.g_lookup = Matrix([[5, 6, 7, 8], [5, 6, 7, 8], [5, 6, 7, 8], [5, 6, 7, 8], [5, 6, 7, 8]])

mca2 = from_workspace(ws.mc_antenna)
mca2.set(mca)

assert mca2 == mca
