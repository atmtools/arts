from pyarts.classes.Rational import Rational


# Get a workspace
r = Rational()

r.set(1.5)
assert r == r
print(r)
r.set(r + 3)
assert r == r
print(r)
r.set(r/2)
assert r == r
print(r)
r.set(r-Rational(5, 7))
assert r == r
print(r)
r.set(r**2)
assert r == r
print(r)
r.set(r**2.5)
assert r == r
print(r)
r.set(3/2)
assert r == r
print(r)

r2 = Rational()
r.savexml("tmp.r.xml", "binary")
r2.readxml("tmp.r.xml")
assert r == r2
