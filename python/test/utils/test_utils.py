# -*- coding: utf-8 -*-
"""Testing the functions in arts.utils.
"""
import warnings
import numpy
from time import sleep
from datetime import timedelta

import pytest

from pyarts import utils


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
