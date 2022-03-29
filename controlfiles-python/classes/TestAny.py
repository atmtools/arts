import pyarts.pyarts_cpp as cxx

try:
    x = cxx.Any()
    assert False
except TypeError:
    pass
except:
    assert False, "Any should do nothing, and not be initialized by anything"

