import pyarts.pyarts_cpp as cxx
import test_functions as test

x = cxx.String("ho")
test.io(x, delete=True)

assert "ho" == x
assert "h" == x[0]

x[1] = 'i'
assert "hi" == x
assert hash("hi") == hash(x)


