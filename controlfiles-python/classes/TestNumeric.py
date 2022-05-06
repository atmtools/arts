import pyarts.arts as cxx
import test_functions as test

x = cxx.Numeric(0)
test.io(x, delete=True)

assert x.val == 0
x.val = x + 2
assert x.val == 2
