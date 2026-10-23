"""
Time is read as RFC 3339, with a space or "T" between the date and the time of day,
where times without an offset from UTC are UTC.
"""

import datetime as dt

import pyarts3 as pyarts

Time = pyarts.arts.Time

same = [
    "2000-01-01 12:00:00",
    "2000-01-01T12:00:00",
    "2000-01-01t12:00:00",
    "2000-01-01 12:00:00Z",
    "2000-01-01T12:00:00z",
    "2000-01-01 12:00:00+00:00",
    "2000-01-01T21:00:00+09:00",
    "2000-01-01T07:30:00-04:30",
]
for s in same:
    assert Time(s) == Time(same[0]), f"{s!r} is not the same time as {same[0]!r}"

assert Time("2000-01-01T21:00:00+09:00").utc_offset == dt.timedelta(hours=9)

d = dt.datetime(1850, 6, 1, 12, 34, 56, 789000, tzinfo=dt.timezone(dt.timedelta(hours=-3)))
assert Time(d.isoformat()) == Time(d), f"{d.isoformat()!r} is not {d!r}"

for s in same:
    t = Time(s)
    u = Time(str(t))
    assert u == t and u.utc_offset == t.utc_offset, f"{str(t)!r} does not read back as {s!r}"

for bad in [
    "2000-01-01",
    "2000-01-01T12:00",
    "2000-01-01 12:00:00+0900",
    "2000-01-01X12:00:00",
    "2001-02-29T00:00:00",
]:
    try:
        Time(bad)
    except Exception:
        continue
    raise AssertionError(f"{bad!r} must not be read as a time")
