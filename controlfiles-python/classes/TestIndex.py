import pyarts.arts as cxx
import test_functions as test

x = cxx.Index(0)
test.io(x, delete=True)

assert x.val == 0
x.val = x + 3
assert x.val == 3
