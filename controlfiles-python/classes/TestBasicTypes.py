from pyarts.workspace import Workspace
from pyarts.classes.BasicTypes import Numeric, Index, String
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()

# Get the agenda
s = from_workspace(ws.output_file_format)
i = from_workspace(ws.abs_lookup_is_adapted)
n = from_workspace(ws.g0)

assert isinstance(s, String), "Bad read"
assert isinstance(n, Numeric), "Bad read"
assert isinstance(i, Index), "Bad read"

s.val = "Hej"
assert s == "Hej", "Bad write"
s.set(String("Help"))
assert s == "Help", "Bad write"

i.val = 3
assert i == 3, "Bad write"
i.set(5)
assert i == 5, "Bad write"

n.val = 3.5
assert n == 3.5, "Bad write"
n.set(3.14)
assert n == 3.14, "Bad write"

s2 = String(0)
s.savexml("tmp.s.xml", "binary")
s2.readxml("tmp.s.xml")
assert s == s2

i2 = Index(0)
i.savexml("tmp.i.xml", "binary")
i2.readxml("tmp.i.xml")
assert i == i2

n2 = Numeric(0)
n.savexml("tmp.n.xml", "binary")
n2.readxml("tmp.n.xml")
assert n == n2

