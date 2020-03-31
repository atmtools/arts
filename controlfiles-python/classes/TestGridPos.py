from pyarts.classes.GridPos import GridPos

# Get the agenda
gp = GridPos()

assert not gp
gp.fd = 0.3, 0.7
assert gp.fd

gp.idx = -1
assert not gp
