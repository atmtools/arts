from pyarts.classes.RadiationVector import RadiationVector
from pyarts.classes.TransmissionMatrix import TransmissionMatrix
from pyarts.classes.Matrix import Matrix
from pyarts.classes.Vector import Vector

tm = TransmissionMatrix(4, 4)
tm[0] = Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
tm[1] = Matrix([[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
tm[2] = Matrix([[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,1]])
tm[3] = Matrix([[1,0,0,0],[0,1,1,0],[0,1,1,0],[0,0,0,0]])

tm2 = TransmissionMatrix()
tm2.set(tm)
assert tm == tm2

rv = RadiationVector(4, 4)
rv[0] = Vector([0,0,0,0])
rv[1] = Matrix([1,0,0,0])
rv[2] = Matrix([0,1,0,0])
rv[3] = Matrix([0,0,1,1])

rv2 = RadiationVector()
rv2.set(rv)
assert rv == rv2

tm3= TransmissionMatrix()
rv3 = RadiationVector()

tm.savexml("tmp.tm.xml", "binary")
rv.savexml("tmp.rv.xml", "binary")

tm3.readxml("tmp.tm.xml")
rv3.readxml("tmp.rv.xml")

assert tm == tm3
assert rv == rv3
