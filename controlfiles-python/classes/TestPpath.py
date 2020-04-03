from pyarts.workspace import Workspace
from pyarts.classes.Ppath import Ppath
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()
ws.atmosphere_dim = 1
ws.cloudbox_on = 0
ws.ppath_inside_cloudbox_do = 0
cloudbox_limits = from_workspace(ws.cloudbox_limits)
z_field = from_workspace(ws.z_field)
z_field.data = [0, 1, 2, 3]
z_field.data = z_field.data.reshape(4,1,1)
z_surface = from_workspace(ws.z_surface)
z_surface.data = [0]
rte_pos = from_workspace(ws.rte_pos)
rte_pos.data = [4]
rte_los = from_workspace(ws.rte_los)
rte_los.data = [91]
ws.ppath_lmax = 0.1
ws.ppathPlaneParallel()

# Get the path
path = from_workspace(ws.ppath)
path2 = Ppath()
path2.set(path)

assert isinstance(path, Ppath), "Bad read"
assert path, "Bad read"
assert path2 == path

path3 = Ppath()
path.savexml("tmp.pp.xml", "binary")
path3.readxml("tmp.pp.xml")
assert path == path3
