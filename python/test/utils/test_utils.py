# -*- coding: utf-8 -*-
"""Testing the functions in arts.utils.
"""
import warnings
import numpy
import xarray
from time import sleep
from datetime import timedelta

import pytest

from arts import utils


class TestUtils:
    """Testing the arts.utils functions."""
    def test_deprecated(self):
        """Test deprecation warning."""
        @utils.deprecated
        def func():
            pass

        with warnings.catch_warnings():
            warnings.simplefilter('error')
            with pytest.raises(DeprecationWarning):
                func()

    def test_image2mpeg(self):
        """Test the behavior when no files are found."""
        with pytest.raises(Exception):
            utils.image2mpeg(glob='', outfile='foo.mp4')

    def test_unique(self):
        assert utils.unique([0, 5, 1, 2, 0, 3, 1]) == [0, 5, 1, 2, 3]

    def test_undo_xarray_floatification(self):
        ds = xarray.Dataset(
            {"a": (["x"], numpy.array([1, 2, 3], dtype="f4")),
             "b": (["x"], numpy.array([2.0, 3.0, 4.0])),
             "c": (["x"], numpy.array(["2010-01-01", "2010-01-02",
                                       "2010-01-03"], dtype="M8"))})
        ds["a"].encoding = {"dtype": numpy.dtype("i4"),
                            "_FillValue": 1234}
        # c should NOT be converted because it's a time
        ds["c"].encoding = {"dtype": numpy.dtype("i8"),
                            "_FillValue": 12345}
        ds2 = utils.undo_xarray_floatification(ds)
        assert ds is not ds2  # has to be a copy
        assert ds["a"].encoding == ds2["a"].encoding
        assert numpy.allclose(ds["a"], ds2["a"])
        assert ds2["a"].dtype == ds2["a"].encoding["dtype"]
        assert (ds2["c"] == ds["c"]).all()
        assert ds2["c"].dtype == ds["c"].dtype
        assert ds2["b"].dtype == ds["b"].dtype

    def test_Timer_methods(self):
        """Test initialization and basic methods of the `Timer` class."""
        t = utils.Timer(verbose=False)
        t.start()
        sleep(0.05)
        dt = t.stop()

        assert timedelta(seconds=0.03) < dt < timedelta(seconds=0.07)

    def test_Timer_block(self):
        """Test usage of a `Timer` object in a with-block."""
        with utils.Timer(verbose=False):
            ...

    def test_Timer_decorator(self):
        """Test usage of a `Timer` object as function decorator."""
        @utils.Timer(verbose=False)
        def foo():
            return

        foo()
