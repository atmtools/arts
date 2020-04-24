from math import isclose

from pyarts.workspace import Workspace
from pyarts.classes.Time import Time
from pyarts.classes import from_workspace

ws = Workspace()
ws.create_variable("Time", "start")
ws.create_variable("Time", "end")
ws.create_variable("Time", "time")
ws.create_variable("Numeric", "dt")
start_time = from_workspace(ws.start)
end_time = from_workspace(ws.end)
cur_time = from_workspace(ws.time)
dt = from_workspace(ws.dt)

# Test Now and Sleep
ws.Now(ws.start)
ws.Sleep(1)
ws.Now(ws.end)
ws.Duration(ws.dt, ws.start, ws.end)
assert isclose(dt, 1, abs_tol=1e-3), \
    "Slept for one second but time differs by more than 1 ms"

# Test SleepUntil
ws.Now(ws.start)
end_time.sec = start_time.sec + 1
ws.SleepUntil(ws.end)
ws.Now(ws.time)
assert isclose(cur_time.sec - start_time.sec, 1, abs_tol=1e-3), \
    "Slept for what should be one second but is off by more than 1 ms..."

# Test set and equality
t = Time()
t.set(start_time)
assert t == start_time, "Times set to one-another do not agree"

# Test IO
t.savexml("t.tmp.xml")
t2 = Time()
t2.readxml("t.tmp.xml")
assert isclose(t.sec, t2.sec), "Time diverges after IO"
