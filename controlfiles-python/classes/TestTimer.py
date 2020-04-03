from time import sleep
from math import isclose

from pyarts.workspace import Workspace
from pyarts.classes.Timer import Timer
from pyarts.classes import from_workspace

ws = Workspace()
ws.timerStart()

time = 1
sleep(time)

ws.timerStop()

t = from_workspace(ws.timer)
assert isinstance(t, Timer)

if t.supported:
    # Test that we are within 2.5 ticks of the true answer
    assert isclose((t.realtime_end - t.realtime_start) / t.tick, time,
                   rel_tol = 2.5 / t.tick)

# t2 = Timer()
# t.savexml("tmp.t.xml", "binary")
# t2.readxml("tmp.t.xml")
# assert t == t2
