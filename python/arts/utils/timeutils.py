# -*- coding: utf-8 -*-

"""This module contains functions for python's datetime/timedelta objects
"""
import functools
import time
from datetime import datetime, timedelta
from numbers import Number

import netCDF4
import numpy as np
import pandas as pd

__all__ = [
    "date2num",
    "num2date",
    "set_time_resolution",
    "to_datetime",
    "to_timedelta",
    "Timer",
]


def set_time_resolution(datetime_obj, resolution):
    """Set the resolution of a python datetime object.

    Args:
        datetime_obj: A python datetime object.
        resolution: A string indicating the required resolution.

    Returns:
        A datetime object truncated to *resolution*.

    Examples:

    .. code-block:: python

        from arts.utils.time import set_time_resolution, to_datetime

        dt = to_datetime("2017-12-04 12:00:00")
        # datetime.datetime(2017, 12, 4, 12, 0)

        new_dt = set_time_resolution(dt, "day")
        # datetime.datetime(2017, 12, 4, 0, 0)

        new_dt = set_time_resolution(dt, "month")
        # datetime.datetime(2017, 12, 1, 0, 0)
    """
    if resolution == "year":
        return set_time_resolution(datetime_obj, "day").replace(month=1, day=1)
    elif resolution == "month":
        return set_time_resolution(datetime_obj, "day").replace(day=1)
    elif resolution == "day":
        return datetime_obj.replace(hour=0, minute=0, second=0, microsecond=0)
    elif resolution == "hour":
        return datetime_obj.replace(minute=0, second=0, microsecond=0)
    elif resolution == "minute":
        return datetime_obj.replace(second=0, microsecond=0)
    elif resolution == "second":
        return datetime_obj.replace(microsecond=0)
    elif resolution == "millisecond":
        return datetime_obj.replace(
            microsecond=int(datetime_obj.microsecond / 1000) * 1000
        )
    else:
        raise ValueError("Cannot set resolution to '%s'!" % resolution)


def to_datetime(obj):
    """Convert an object to a python datetime object.

    Args:
        obj: Can be a string with time information, a numpy.datetime64 or a
            pandas.Timestamp object.

    Returns:
        A python datetime object.

    Examples:

    .. code-block:: python

        dt = to_datetime("2017-12-04 12:00:00")
        # dt is datetime.datetime(2017, 12, 4, 12, 0)
    """

    if isinstance(obj, datetime):
        return obj
    else:
        return pd.to_datetime(obj).to_pydatetime()


def to_timedelta(obj, numbers_as=None):
    """Convert an object to a python datetime object.

    Args:
        obj: Can be a string with time information, a number, a
            numpy.timedelta64 or a pandas.Timedelta object.
        numbers_as: A string that indicates how numbers should be
            interpreted. Allowed values are *weeks*, *days*, *hours*,
            *minutes*, *seconds*, *milliseconds* and *microseconds.

    Returns:
        A python datetime object.

    Examples:

    .. code-block:: python

        # A timedelta object with 200 seconds
        t = to_timedelta("200 seconds")

        # A timedelta object with 24 days
        t = to_timedelta(24, numbers_as="days")
    """

    if numbers_as is None:
        numbers_as = "seconds"

    if isinstance(obj, timedelta):
        return obj
    elif isinstance(obj, Number):
        return timedelta(**{numbers_as: int(obj)})
    else:
        return pd.to_timedelta(obj).to_pytimedelta()


unit_mapper = {
    "nanoseconds": "ns",
    "microseconds": "us",
    "milliseconds": "ms",
    "seconds": "s",
    "hours": "h",
    "minutes": "m",
    "days": "D",
}


class InvalidUnitString(Exception):
    def __init__(self, *args, **kwargs):
        super(InvalidUnitString, self).__init__(*args, **kwargs)


def date2num(dates, units, calendar=None):
    """Convert an array of integer into datetime objects.

    This function optimizes the date2num function of python-netCDF4 if the
    standard calendar is used.

    Args:
        dates: Either an array of numpy.datetime64 objects (if standard
            gregorian calendar is used), otherwise an array of python
            datetime objects.
        units: A string with the format "{unit} since {epoch}",
            e.g. "seconds since 1970-01-01T00:00:00".
        calendar: (optional) Standard is gregorian. If others are used,
            netCDF4.num2date will be called.

    Returns:
        An array of integers.
    """
    if calendar is None:
        calendar = "gregorian"
    else:
        calendar = calendar.lower()

    if calendar != "gregorian":
        return netCDF4.date2num(dates, units, calendar)

    try:
        unit, epoch = units.split(" since ")
    except ValueError:
        raise InvalidUnitString("Could not convert to numeric values!")

    converted_data = \
        dates.astype("M8[%s]" % unit_mapper[unit]).astype("int")

    # numpy.datetime64 cannot read certain time formats while pandas can.
    epoch = pd.Timestamp(epoch).to_datetime64()

    if epoch != np.datetime64("1970-01-01"):
        converted_data -= np.datetime64("1970-01-01") - epoch
    return converted_data


def num2date(times, units, calendar=None):
    """Convert an array of integers into datetime objects.

    This function optimizes the num2date function of python-netCDF4 if the
    standard calendar is used.

    Args:
        times: An array of integers representing timestamps.
        units: A string with the format "{unit} since {epoch}",
            e.g. "seconds since 1970-01-01T00:00:00".
        calendar: (optional) Standard is gregorian. If others are used,
            netCDF4.num2date will be called.

    Returns:
        Either an array of numpy.datetime64 objects (if standard gregorian
        calendar is used), otherwise an array of python datetime objects.
    """
    try:
        unit, epoch = units.split(" since ")
    except ValueError:
        raise InvalidUnitString("Could not convert to datetimes!")

    if calendar is None:
        calendar = "gregorian"
    else:
        calendar = calendar.lower()

    if calendar != "gregorian":
        return netCDF4.num2date(times, units, calendar).astype(
            "M8[%s]" % unit_mapper[unit])

    # Numpy uses the epoch 1970-01-01 natively.
    converted_data = times.astype("M8[%s]" % unit_mapper[unit])

    # numpy.datetime64 cannot read certain time formats while pandas can.
    epoch = pd.Timestamp(epoch).to_datetime64()

    # Maybe there is another epoch used?
    if epoch != np.datetime64("1970-01-01"):
        converted_data -= np.datetime64("1970-01-01") - epoch
    return converted_data


class Timer:
    """Provide a simple time profiling utility

    Parameters:
        verbose (bool): Print measured duration after stopping the timer.
        info (str): Allows to add additional information to output.
            The given string is printed before the measured time.
            If `None`, default information is added depending on the use case.

    Returns:
        datetime.timedelta: The duration between start and end time.

    Examples:
        Timer in with statement:

        >>> import time
        >>> with Timer():
        ...     time.sleep(1)
        elapsed time: 0:00:01.003186

        Timer as object (allows to store :class:`datetime.timedelta`):

        >>> import time
        >>> t = Timer().start()
        >>> time.sleep(1)
        >>> dt = t.stop()
        elapsed time: 0:00:01.004756

        As function decorator:

        >>> @Timer()
        ... def own_function(s):
        ...     import time
        ...     time.sleep(s)
        >>> own_function(1)
        own_function: 0:00:01.004667

        Use it in format strings:

        >>> from arts.utils import Timer
        >>> timer = Timer().start()
        >>> print(f"{timer} elapsed")
        0:00:00.000111 hours elapsed
    """
    def __init__(self, info=None, verbose=True):
        """Create a timer object."""
        self.verbose = verbose
        self.info = info

        self.starttime = None
        self.endtime = None

    def __call__(self, func):
        """Allows to use a Timer object as a decorator."""
        # When no additional information is given, add the function name if
        # Timer is used as decorator.
        if self.info is None:
            self.info = func.__name__

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            with self:
                # Call the original function in a Timer context.
                return func(*args, **kwargs)

        return wrapper

    def __enter__(self):
        return self.start()

    def __exit__(self, *args):
        self.stop()

    def __repr__(self):
        return repr(self.elapsed)

    def __str__(self):
        return str(self.elapsed) + " hours"

    @property
    def elapsed(self):
        """Get the elapsed time as timedelta object"""
        if self.starttime is None:
            raise ValueError("Timer has not been started yet!")

        return timedelta(seconds=time.time() - self.starttime)

    def start(self):
        """Start timer."""
        self.starttime = time.time()
        return self

    def stop(self):
        """Stop timer and print info message

        The info message will be only printed if `Timer.verbose` is *true*.

        Returns:
            A timedelta object.
        """
        if self.starttime is None:
            raise ValueError("Timer has not been started yet!")

        self.endtime = time.time()

        dt = timedelta(seconds=self.endtime - self.starttime)

        # If no additional information is specified add default information
        # to make the output more readable.
        if self.info is None:
            self.info = 'elapsed time'

        if self.verbose:
            # Connect additional information and measured time for output.
            print('{info}: {time}'.format(info=self.info, time=dt))

        return dt
