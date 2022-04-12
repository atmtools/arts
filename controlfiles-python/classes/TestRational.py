import pyarts.pyarts_cpp as cxx
import test_functions as test

x = cxx.Rational()
test.io(x, delete=True)

y = x + 1
x += 1
assert y == x

y = x * 2
x *= 2
assert y == x

y = x / 3
x /= 3
assert y == x

y = x - 4
x -= 4
assert y == x

