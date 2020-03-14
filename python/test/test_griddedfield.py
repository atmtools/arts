# -*- encoding: utf-8 -*-
import os
from tempfile import mkstemp

import numpy as np
import pytest

from pyarts import griddedfield, xml


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
        1: griddedfield.GriddedField1,
        2: griddedfield.GriddedField2,
        3: griddedfield.GriddedField3,
        4: griddedfield.GriddedField4,
        5: griddedfield.GriddedField5,
        6: griddedfield.GriddedField6,
    }.get(dim)


class TestGriddedFieldUsage:
    ref_dir = os.path.join(os.path.dirname(__file__), "reference", "")

    @pytest.mark.parametrize("dim", range(1, 7))
    def test_check_init(self, dim):
        """Test initialisation of GriddedFields."""
        cls = _get_griddedfield_type(dim)
        assert cls().dimension == dim

    def test_check_dimension1(self):
        """Test if grid and data dimension agree (positive)."""
        gf3 = griddedfield.GriddedField3()
        gf3.grids = [np.arange(5), np.arange(5), []]
        gf3.gridnames = ["A", "B", "C"]
        gf3.data = np.ones((5, 5, 1))
        assert gf3.check_dimension() is True

    def test_check_dimension2(self):
        """Test if grid and data dimension agree (negative)."""
        gf3 = griddedfield.GriddedField3()
        gf3.grids = [np.arange(5), np.arange(5), []]
        gf3.gridnames = ["A", "B", "C"]

        # It shold not be allowed to set a Matrix as data in a GriddedField3.
        with pytest.raises(TypeError):
            gf3.data = np.ones((5, 5))

    def test_data(self):
        """Test setting and getting of data. """
        reference = np.random.randn(10, 10, 10)
        gf3 = griddedfield.GriddedField3()
        gf3.data = reference
        assert np.array_equal(gf3.data, reference)

    def test_name_setter(self):
        """Test name setter and getter."""
        reference = 'TestName'
        gf = griddedfield.GriddedField1()
        gf.name = reference
        assert gf.name == reference

    def _setup_gf2(self):
        """Helper for test_to_dict and test_to_xarray"""
        gf2 = griddedfield.GriddedField2()
        gf2.grids = [np.ones(5), np.zeros(5)]
        gf2.gridnames = ["ones", "zeros"]
        gf2.data = np.ones((5, 5))
        gf2.name = "semprini"
        return gf2

    def test_to_dict(self):
        """Test the conversion into a dictionary."""
        gf2 = self._setup_gf2()
        d = gf2.to_dict()

        res = (np.array_equal(d['ones'], np.ones(5)) and
               np.array_equal(d['zeros'], np.zeros(5)) and
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
            gf = griddedfield.GriddedField1()
            gf.name = nametype

    def test_shape(self):
        """Test return of data shape."""
        gf3 = xml.load(self.ref_dir + 'GriddedField3.xml')

        assert gf3.shape == gf3.data.shape == (2, 2, 2)

    def test_data_subscription(self):
        """Test direct data subscription."""
        gf3 = xml.load(self.ref_dir + 'GriddedField3.xml')

        assert gf3[0, 1, 0] == gf3.data[0, 1, 0]

    def test_slicing(self):
        """Test GriddedField slicing."""
        gf3 = xml.load(self.ref_dir + 'GriddedField3.xml')

        # Create new GriddedField which is a sliced subset of the initial one.
        gf_sliced = gf3.extract_slice(slice(1, None), axis=1)

        assert np.allclose(gf3.data[:, 1:, :], gf_sliced.data)

    def test_repr(self):
        """Test string represenation of GriddedField objects."""
        str(xml.load(self.ref_dir + 'GriddedField3.xml'))

    def test_repr_empty(self):
        """Test string represenation of empty GriddedField objects."""
        str(griddedfield.GriddedField1())

    def test_get(self):
        """Test the get method for named fields."""
        gf1 = griddedfield.GriddedField1(
            grids=[['foo', 'bar']],
            data=np.array([42, 13]),
        )

        assert gf1.get('foo') == np.array([42])

    def test_get_default(self):
        """Test the GriddedField.get() behavior for non-existing fieldnames."""
        gf1 = griddedfield.GriddedField1(
            grids=[['dummy']],
            data=np.array([0]),
        )

        # Return given default, if a name is not existing.
        assert gf1.get('nonexisting', 42) == 42

        # If no default is specified, return `None`.
        assert gf1.get('nonexisting') is None

    def test_get_keepdims(self):
        """Test the dimension handling of the GriddedField.get()."""
        gf1 = griddedfield.GriddedField1(
            grids=[['foo', 'bar']],
            data=np.array([42, 13]),
        )

        assert gf1.get('foo').shape == (1,)
        assert gf1.get('foo', keep_dims=False).shape == tuple()

    def test_get_nofieldnames(self):
        """Test behavior if first grids is not ArrayOfString."""
        gf1 = griddedfield.GriddedField1(
            grids=[[0]],
            data=np.array([0]),
        )

        with pytest.raises(TypeError):
            gf1.get(0)

    def test_scaling(self):
        """Test the scaling of data in named fields."""
        gf1 = griddedfield.GriddedField1(
            grids=[['first_field', 'second_field']],
            data=np.array([1., 1.]),
        )

        gf1.scale('second_field', 0.1)

        # Check if values if *only* values of the second fields are scaled.
        assert gf1.data[0] == np.array([1])
        assert gf1.data[1] == np.array([0.1])

    def test_integer_scaling(self):
        """Test the scaling of integer data in named fields."""
        gf1 = griddedfield.GriddedField1(
            grids=[['first_field', 'second_field']],
            data=np.array([1, 1]),
        )

        gf1.scale('second_field', 0.1)

        # Check if values if *only* values of the second fields are scaled.
        assert gf1.data[0] == np.array([1])
        assert gf1.data[1] == np.array([0.1])

    def test_set(self):
        """Test the set method for named fields."""
        gf1 = griddedfield.GriddedField1(
            grids=[['zero', 'one']],
            data=np.array([0, 0]),
        )

        gf1.set('one', 1)

        assert gf1.data[1] == np.array([1])


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
        assert np.array_equal(test_data, reference)

    def test_load_dimension(self):
        """Load reference XML file for GriddedField3 and run check."""
        gf3 = xml.load(self.ref_dir + 'GriddedField3.xml')
        assert gf3.check_dimension()

    def test_equality(self):
        """Check the equality of two GriddedField objects."""
        # Create two different objects with same content.
        a = xml.load(self.ref_dir + 'GriddedField3.xml')
        b = xml.load(self.ref_dir + 'GriddedField3.xml')

        assert a == b

    def test_equality_empty(self):
        """Check the equality of two empty GriddedField objects."""
        # Create two different objects with same content.
        a = griddedfield.GriddedField3()
        b = griddedfield.GriddedField3()

        assert a == b

    def test_nonequality(self):
        """Check the non-equality of two GriddedField objects."""
        # Create two different objects with same content.
        a = xml.load(self.ref_dir + 'GriddedField3.xml')
        b = xml.load(self.ref_dir + 'GriddedField3.xml')

        a.name = 'foo'
        b.name = 'bar'

        assert a != b

    def test_copy(self):
        """Test copying of GriddedFields."""
        a = xml.load(self.ref_dir + 'GriddedField3.xml')
        b = a.copy()

        # GriddedFields should be equal but not the same object.
        assert a is not b and a == b

    def test_deepcopy(self):
        """Test deepcopying of GriddedField attributes."""
        a = xml.load(self.ref_dir + 'GriddedField3.xml')
        b = a.copy()

        # Grids should not be the same object.
        assert a.grids is not b.grids

    def test_from_xarray(self):
        a = xml.load(self.ref_dir + 'GriddedField3.xml')
        a.dataname = 'Testdata'
        da = a.to_xarray()
        b = griddedfield.GriddedField3.from_xarray(da)
        assert a == b


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
