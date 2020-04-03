import numpy as np

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

n = int(np.random.permutation(np.linspace(1,10,10))[0])
t1.data = np.random.normal(size=(n,))
t2.data = np.random.normal(size=(n, n))
t3.data = np.random.normal(size=(n, n, n))
t4.data = np.random.normal(size=(n, n, n, n))
t5.data = np.random.normal(size=(n, n, n, n, n))
t6.data = np.random.normal(size=(n, n, n, n, n, n))
t7.data = np.random.normal(size=(n, n, n, n, n, n, n))

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

assert isinstance(t1 @ t1, float)
assert isinstance(t2 @ t2, np.ndarray)
assert isinstance(t3 @ t3, np.ndarray)
assert isinstance(t4 @ t4, np.ndarray)
assert isinstance(t5 @ t5, np.ndarray)
assert isinstance(t6 @ t6, np.ndarray)
assert isinstance(t7 @ t7, np.ndarray)

assert isinstance(t1 * t1, Vector)
assert isinstance(t1 / t1, Vector)
assert isinstance(t1 - t1, Vector)
assert isinstance(t1 + t1, Vector)
assert isinstance(t2 * t1, Matrix)
assert isinstance(t2 / t1, Matrix)
assert isinstance(t2 - t1, Matrix)
assert isinstance(t2 + t1, Matrix)
assert isinstance(t3 * t1, Tensor3)
assert isinstance(t3 / t1, Tensor3)
assert isinstance(t3 - t1, Tensor3)
assert isinstance(t3 + t1, Tensor3)
assert isinstance(t4 * t1, Tensor4)
assert isinstance(t4 / t1, Tensor4)
assert isinstance(t4 - t1, Tensor4)
assert isinstance(t4 + t1, Tensor4)
assert isinstance(t5 * t1, Tensor5)
assert isinstance(t5 / t1, Tensor5)
assert isinstance(t5 - t1, Tensor5)
assert isinstance(t5 + t1, Tensor5)
assert isinstance(t6 * t1, Tensor6)
assert isinstance(t6 / t1, Tensor6)
assert isinstance(t6 - t1, Tensor6)
assert isinstance(t6 + t1, Tensor6)
assert isinstance(t7 * t1, Tensor7)
assert isinstance(t7 / t1, Tensor7)
assert isinstance(t7 - t1, Tensor7)
assert isinstance(t7 + t1, Tensor7)

assert isinstance(t1 * 2, Vector)
assert isinstance(t1 / 2, Vector)
assert isinstance(t1 - 2, Vector)
assert isinstance(t1 + 2, Vector)
assert isinstance(t2 * 2, Matrix)
assert isinstance(t2 / 2, Matrix)
assert isinstance(t2 - 2, Matrix)
assert isinstance(t2 + 2, Matrix)
assert isinstance(t3 * 2, Tensor3)
assert isinstance(t3 / 2, Tensor3)
assert isinstance(t3 - 2, Tensor3)
assert isinstance(t3 + 2, Tensor3)
assert isinstance(t4 * 2, Tensor4)
assert isinstance(t4 / 2, Tensor4)
assert isinstance(t4 - 2, Tensor4)
assert isinstance(t4 + 2, Tensor4)
assert isinstance(t5 * 2, Tensor5)
assert isinstance(t5 / 2, Tensor5)
assert isinstance(t5 - 2, Tensor5)
assert isinstance(t5 + 2, Tensor5)
assert isinstance(t6 * 2, Tensor6)
assert isinstance(t6 / 2, Tensor6)
assert isinstance(t6 - 2, Tensor6)
assert isinstance(t6 + 2, Tensor6)
assert isinstance(t7 * 2, Tensor7)
assert isinstance(t7 / 2, Tensor7)
assert isinstance(t7 - 2, Tensor7)
assert isinstance(t7 + 2, Tensor7)

assert isinstance(2 * t1, Vector)
assert isinstance(2 / t1, Vector)
assert isinstance(2 - t1, Vector)
assert isinstance(2 + t1, Vector)
assert isinstance(2 * t2, Matrix)
assert isinstance(2 / t2, Matrix)
assert isinstance(2 - t2, Matrix)
assert isinstance(2 + t2, Matrix)
assert isinstance(2 * t3, Tensor3)
assert isinstance(2 / t3, Tensor3)
assert isinstance(2 - t3, Tensor3)
assert isinstance(2 + t3, Tensor3)
assert isinstance(2 * t4, Tensor4)
assert isinstance(2 / t4, Tensor4)
assert isinstance(2 - t4, Tensor4)
assert isinstance(2 + t4, Tensor4)
assert isinstance(2 * t5, Tensor5)
assert isinstance(2 / t5, Tensor5)
assert isinstance(2 - t5, Tensor5)
assert isinstance(2 + t5, Tensor5)
assert isinstance(2 * t6, Tensor6)
assert isinstance(2 / t6, Tensor6)
assert isinstance(2 - t6, Tensor6)
assert isinstance(2 + t6, Tensor6)
assert isinstance(2 * t7, Tensor7)
assert isinstance(2 / t7, Tensor7)
assert isinstance(2 - t7, Tensor7)
assert isinstance(2 + t7, Tensor7)

assert isinstance(t2 @ t1, np.ndarray)
assert isinstance(t3 @ t2, np.ndarray)
assert isinstance(t4 @ t3, np.ndarray)
assert isinstance(t5 @ t4, np.ndarray)
assert isinstance(t6 @ t5, np.ndarray)
assert isinstance(t7 @ t6, np.ndarray)

t7 @= t6
assert isinstance(t7, Tensor7)
t6 @= t5
assert isinstance(t6, Tensor6)
t5 @= t4
assert isinstance(t5, Tensor5)
t4 @= t3
assert isinstance(t4, Tensor4)
t3 @= t2
assert isinstance(t3, Tensor3)
t2 @= t1
assert isinstance(t2, Matrix)
t1 @= t1
assert isinstance(t1, Vector)

t7.data = t1
assert (t1 == t7).all()
t7.data = t2
assert (t2 == t7).all()
t7.data = t3
assert (t3 == t7).all()
t7.data = t4
assert (t4 == t7).all()
t7.data = t5
assert (t5 == t7).all()
t7.data = t6
assert (t6 == t7).all()

t12 = Vector()
t22 = Matrix()
t32 = Tensor3()
t42 = Tensor4()
t52 = Tensor5()
t62 = Tensor6()
t72 = Tensor7()

t12.set(t1)
t22.set(t2)
t32.set(t3)
t42.set(t4)
t52.set(t5)
t62.set(t6)
t72.set(t7)

assert t12 == t1
assert t22 == t2
assert t32 == t3
assert t42 == t4
assert t52 == t5
assert t62 == t6
assert t72 == t7

t13 = Vector()
t1.savexml("tmp.mt.xml", "binary")
t13.readxml("tmp.mt.xml")
assert t1 == t13

t23 = Matrix()
t2.savexml("tmp.mt.xml", "binary")
t23.readxml("tmp.mt.xml")
assert t2 == t23

t33 = Tensor3()
t3.savexml("tmp.mt.xml", "binary")
t33.readxml("tmp.mt.xml")
assert t3 == t33

t43 = Tensor4()
t4.savexml("tmp.mt.xml", "binary")
t43.readxml("tmp.mt.xml")
assert t4 == t43

t53 = Tensor5()
t5.savexml("tmp.mt.xml", "binary")
t53.readxml("tmp.mt.xml")
assert t5 == t53

t63 = Tensor6()
t6.savexml("tmp.mt.xml", "binary")
t63.readxml("tmp.mt.xml")
assert t6 == t63

t73 = Tensor7()
t7.savexml("tmp.mt.xml", "binary")
t73.readxml("tmp.mt.xml")
assert t7 == t73
