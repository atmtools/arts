from pyarts.classes.GridPos import GridPos

# Get the agenda
gp = GridPos()

assert gp
gp.fd = 0.3, 0.7
assert gp.fd

gp.idx = -1
assert not gp

gp2 = GridPos()
gp2.set(gp)

assert gp2 == gp
