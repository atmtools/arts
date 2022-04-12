from pyarts.classes.GridPos import GridPos

# Get the agenda
gp = GridPos()

assert gp
gp.fd = 0.3, 0.7
assert gp.fd

gp.idx = -1
assert not gp
gp.idx = 0
assert gp

gp2 = GridPos()
gp2.set(gp)

assert gp2 == gp

gp3 = GridPos()
gp.savexml("tmp.gp.xml", "binary")
gp3.readxml("tmp.gp.xml")
assert gp == gp3
