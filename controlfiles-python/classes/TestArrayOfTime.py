import pyarts.pyarts_cpp as cxx
import test_functions as test

import datetime as datetime

x = cxx.ArrayOfTime(1, "2017-01-01 15:30:20")
test.io(x, delete=True)
test.array(x)

x = cxx.ArrayOfTime([f"2020-01-{x} 00:00:00" for x in range(10, 30)])

assert x.as_datetime[-1] == datetime.datetime(2020, 1, 29, 0, 0)
