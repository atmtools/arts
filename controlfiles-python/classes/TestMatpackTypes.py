from pyarts.classes.Vector import Vector
from pyarts.classes.Matrix import Matrix
from pyarts.classes.Tensor3 import Tensor3
from pyarts.classes.Tensor4 import Tensor4
from pyarts.classes.Tensor5 import Tensor5
from pyarts.classes.Tensor6 import Tensor6
from pyarts.classes.Tensor7 import Tensor7

t1 = Vector()
t2 = Matrix()
t3 = Tensor3()
t4 = Tensor4()
t5 = Tensor5()
t6 = Tensor6()
t7 = Tensor7()

assert not t1
assert not t2
assert not t3
assert not t4
assert not t5
assert not t6
assert not t7

t1.data = [1, 2, 3]
t2.data = [1, 2, 3]
t3.data = [1, 2, 3]
t4.data = [1, 2, 3]
t5.data = [1, 2, 3]
t6.data = [1, 2, 3]
t7.data = [1, 2, 3]

assert t1
assert t2
assert t3
assert t4
assert t5
assert t6
assert t7

assert t1 == t1
assert t2 == t2
assert t3 == t3
assert t4 == t4
assert t5 == t5
assert t6 == t6
assert t7 == t7

assert (t1.data == [1, 2, 3]).all()
assert (t1.data == t2.data).all()
assert (t1.data == t3.data).all()
assert (t1.data == t4.data).all()
assert (t1.data == t5.data).all()
assert (t1.data == t6.data).all()
assert (t1.data == t7.data).all()

t7.data = t7.data.reshape(3, 1, 1, 1, 1, 1, 1)
assert (t1.data == t7.data.flatten()).all()
