# -*- encoding: utf-8 -*-
import os
from tempfile import mkstemp
from copy import deepcopy as copy

import xarray as xa
import numpy as np
import pytest

from pyarts import xml
import pyarts.arts as cxx


def _create_tensor(n):
    """Create a tensor of dimension n.

    Create a tensor with n dimensions with two entries in each dimension.
    The tensor is filled with increasing integers starting with 0.

    Args:
        n (int): number of dimensions

    Returns:
        np.ndarray: n-dimensional tensor

    """
    return np.arange(2 ** n).reshape(2 * np.ones(n).astype(int))


def _get_griddedfield_type(dim):
    """Return the appropriate `GriddedField` type for given dimension."""
    return {
        1: cxx.GriddedField1,
        2: cxx.GriddedField2,
        3: cxx.GriddedField3,
        4: cxx.GriddedField4,
        5: cxx.GriddedField5,
        6: cxx.GriddedField6,
    }.get(dim)


class TestGriddedFieldUsage:
    ref_dir = os.path.join(os.path.dirname(__file__), "reference", "")

    @pytest.mark.parametrize("dim", range(1, 7))
    def test_check_init(self, dim):
        """Test initialisation of GriddedFields."""
        cls = _get_griddedfield_type(dim)
        assert cls().dim == dim

    def test_check_dimension1(self):
        """Test if grid and data dimension agree (positive)."""
        gf3 = cxx.GriddedField3()
        gf3.grids = [np.arange(5), np.arange(5), []]
        gf3.gridnames = ["A", "B", "C"]
        gf3.data = np.ones((5, 5, 0))
        assert gf3.ok() is True

    def test_data(self):
        """Test setting and getting of data. """
        reference = np.random.randn(10, 10, 10)
        gf3 = cxx.GriddedField3()
        gf3.data = reference
        assert np.array_equal(gf3.data, reference)

    def test_name_setter(self):
        """Test name setter and getter."""
        reference = 'TestName'
        gf = cxx.GriddedField1()
        gf.dataname = reference
        assert gf.dataname == reference

    def _setup_gf2(self):
        """Helper for test_to_dict and test_to_xarray"""
        gf2 = cxx.GriddedField2()
        gf2.grids = [np.ones(5), np.zeros(5)]
        gf2.gridnames = ["ones", "zeros"]
        gf2.data = np.ones((5, 5))
        gf2.dataname = "semprini"
        return gf2

    def test_to_dict(self):
        """Test the conversion into a dictionary."""
        gf2 = self._setup_gf2()
        d = gf2.to_dict()

        res = (np.array_equal(d['coords']['ones']['data'], np.ones(5)) and
               np.array_equal(d['coords']['zeros']['data'], np.zeros(5)) and
               np.array_equal(d['data'], np.ones((5, 5))))

        assert res is True

    def test_to_xarray(self):
        """Test the conversion into xarray DataArray"""
        gf2 = self._setup_gf2()

        da = gf2.to_xarray()

        assert (da.name == "semprini" and
                da.dims == ("ones", "zeros") and
                np.array_equal(da.coords["zeros"], np.zeros(5)) and
                np.array_equal(da.coords["ones"], np.ones(5)) and
                np.array_equal(da.values, np.ones((5, 5))))

    @pytest.mark.parametrize('nametype', [float(), int()])
    def test_name_type(self, nametype):
        """Test if only names of type str are accepted."""
        with pytest.raises(TypeError):
            gf = cxx.GriddedField1()
            gf.dataname = nametype

    def test_shape(self):
        """Test return of data shape."""
        gf3 = xml.load(self.ref_dir + 'GriddedField3.xml')
        assert gf3.shape == gf3.data.shape == (2, 2, 2)

    def test_data_subscription(self):
        """Test direct data subscription."""
        gf3 = xml.load(self.ref_dir + 'GriddedField3.xml')
        assert gf3[0, 1, 0] == gf3.data[0, 1, 0]

    def test_repr(self):
        """Test string represenation of GriddedField objects."""
        str(xml.load(self.ref_dir + 'GriddedField3.xml'))

    def test_repr_empty(self):
        """Test string represenation of empty GriddedField objects."""
        str(cxx.GriddedField1())


class TestGriddedFieldLoad:
    ref_dir = os.path.join(os.path.dirname(__file__), "reference", "")

    def test_load_data(self):
        """Load reference XML file for GriddedField3 and check the data."""
        reference = _create_tensor(3)
        gf3 = xml.load(self.ref_dir + 'GriddedField3.xml')
        test_data = gf3.data
        assert np.array_equal(test_data, reference)

    def test_load_grids(self):
        """Load reference XML file for GriddedField3 and check the grids."""
        reference = [np.arange(2)] * 3
        gf3 = xml.load(self.ref_dir + 'GriddedField3.xml')
        test_data = gf3.grids
        assert all(np.allclose(a, b) for a, b in zip(test_data, reference))

    def test_load_gridnames(self):
        """Load reference XML file for GriddedField3 and check gridnames."""
        reference = ['grid1', 'grid2', 'grid3']
        gf3 = xml.load(self.ref_dir + 'GriddedField3.xml')
        test_data = gf3.gridnames
        
        assert test_data == reference

    def test_load_dimension(self):
        """Load reference XML file for GriddedField3 and run check."""
        gf3 = xml.load(self.ref_dir + 'GriddedField3.xml')
        assert gf3.ok()

    def test_equality(self):
        """Check the equality of two GriddedField objects."""
        # Create two different objects with same content.
        a = xml.load(self.ref_dir + 'GriddedField3.xml')
        b = xml.load(self.ref_dir + 'GriddedField3.xml')

        assert a == b

    def test_equality_empty(self):
        """Check the equality of two empty GriddedField objects."""
        # Create two different objects with same content.
        a = cxx.GriddedField3()
        b = cxx.GriddedField3()

        assert a == b

    def test_nonequality(self):
        """Check the non-equality of two GriddedField objects."""
        # Create two different objects with same content.
        a = xml.load(self.ref_dir + 'GriddedField3.xml')
        b = xml.load(self.ref_dir + 'GriddedField3.xml')

        a.dataname = 'foo'
        b.dataname = 'bar'

        assert a != b

    def test_copy(self):
        """Test copying of GriddedFields."""
        a = xml.load(self.ref_dir + 'GriddedField3.xml')
        b = copy(a)

        # GriddedFields should be equal but not the same object.
        assert a is not b and a == b

    def test_deepcopy(self):
        """Test deepcopying of GriddedField attributes."""
        a = xml.load(self.ref_dir + 'GriddedField3.xml')
        b = copy(a)

        # Grids should not be the same object.
        assert a.grids is not b.grids

    def test_from_xarray(self):
        a = xml.load(self.ref_dir + 'GriddedField3.xml')
        a.dataname = 'Testdata'
        da = a.to_xarray()
        b = cxx.GriddedField3.from_xarray(da)
        assert a == b

        v = xa.DataArray(
            data=np.zeros(shape=(3, 1, 1)),
            dims=["alt", "lat", "lon"],
            coords={"alt": [0., 1, 2], "lat": [0.], "lon": [0.]},
        )
        c = cxx.GriddedField3.from_xarray(v)


class TestGriddedFieldWrite:
    def setup_method(self):
        """Create a temporary file."""
        fd, self.f = mkstemp()
        os.close(fd)

    def teardown_method(self):
        """Delete temporary file."""
        os.remove(self.f)

    @pytest.mark.parametrize("dim", range(1, 7))
    def test_write_load_griddedfield(self, dim):
        gf = _get_griddedfield_type(dim)()
        gf.grids = [np.arange(2)] * dim
        gf.data = _create_tensor(dim)
        xml.save(gf, self.f)

        test_data = xml.load(self.f)

        assert np.array_equal(gf.data, test_data.data)


if __name__ == "__main__":
    x = TestGriddedFieldLoad()
    x.test_from_xarray()
