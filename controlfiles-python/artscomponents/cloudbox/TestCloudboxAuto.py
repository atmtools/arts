from pyarts.workspace import Workspace
from pyarts.arts import Tensor4
import numpy as np

ws = Workspace()
ws.AtmosphereSet1D()
p_grid = np.linspace(10000, 1, 100)
ws.p_grid = p_grid

lowercloudlimit = 5
uppercloudlimit = 20
particle_field = np.zeros([1, 100, 1, 1])
particle_field[0, lowercloudlimit:uppercloudlimit + 1, 0, 0] = 1
particle_field_tensor4 = Tensor4(particle_field)

# Check that default cloudbox_margin=-1 works correctly
ws.cloudboxSetAutomatically(particle_field=particle_field_tensor4)
cloudbox_limits = ws.cloudbox_limits.value
assert cloudbox_limits[0] == 0
assert cloudbox_limits[1] == uppercloudlimit + 1

# Check user-defined cloudbox_margin works correctly
ws.cloudboxSetAutomatically(particle_field=particle_field_tensor4, cloudbox_margin=10)
cloudbox_limits = ws.cloudbox_limits.value
assert cloudbox_limits[0] == lowercloudlimit - 2
assert cloudbox_limits[1] == uppercloudlimit + 1
