# -*- coding: utf-8 -*-
"""Testing the basic ARTS XML functions

This module provides basic functions to test the reading and writing
of ARTS XML files.
"""
import os
from tempfile import mkstemp

import numpy as np
from scipy import sparse
import pytest

from pyarts import xml
from pyarts.catalogues import Sparse


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


def _create_complex_tensor(n):
    """Create a complex tensor of dimension n.

    Create a complex tensor with n dimensions with two entries in each
    dimension. The tensor is filled with increasing integers starting with 0.

    Args:
        n (int): number of dimensions

    Returns:
        np.ndarray: n-dimensional tensor

    """
    return np.arange(2 ** n,
                     dtype=np.complex128).reshape(2 * np.ones(n).astype(int))


def _create_empty_tensor(n):
    """Create an empty tensor of dimension n.

    Create a tensor with n dimensions of size 0.

    Args:
        n (int): number of dimensions

    Returns:
        np.ndarray: n-dimensional tensor

    """
    return np.ndarray((0,) * n)


def _create_sparse(n):
    """Create a n x n diagonal sparse with ones on the diagonal.

    Args:
        n (int): size

    Returns:
        arts.catalogues.Sparse: sparse matrix with ones as diagonals

    """
    a = sparse.diags(np.ones(n), format='csc')
    return Sparse(a)


class TestLoad:
    """Testing the ARTS XML reading functions.

    This class provides functions to test the reading of XML files. For this
    purpose reference files are read and compared to the expexted results.

    Attributes:
        ref_dir (str): absolute path to the reference data directory.

    """
    ref_dir = os.path.join(os.path.dirname(__file__), "reference", "")

    def test_load_index(self):
        """Load reference XML file for ARTS type Index."""
        assert xml.load(self.ref_dir + 'index.xml') == 0

    def test_load_vector(self):
        """Load reference XML file for ARTS type Vector."""
        reference = _create_tensor(1)
        test_data = xml.load(self.ref_dir + 'vector.xml')
        assert np.array_equal(test_data, reference)

    def test_load_vector_binary(self):
        """Load reference binary XML file for ARTS type Vector."""
        reference = _create_tensor(1)
        test_data = xml.load(self.ref_dir + 'vector-bin.xml')
        assert np.array_equal(test_data, reference)

    def test_load_matrix(self):
        """Load reference XML file for ARTS type Matrix."""
        reference = _create_tensor(2)
        test_data = xml.load(self.ref_dir + 'matrix.xml')
        assert np.array_equal(test_data, reference)

    @pytest.mark.parametrize('suffix', ['.xml', '-bin.xml'])
    def test_load_sparse(self, suffix):
        """Load reference XML file for ARTS type Sparse."""
        reference = _create_sparse(10) * 7
        test_data = xml.load(self.ref_dir + 'sparse' + suffix)
        assert np.array_equal(test_data.toarray(), reference.toarray())

    @pytest.mark.parametrize('n', range(3, 8))
    def test_load_tensor(self, n):
        """Load tensor of dimension n and compare data to reference.

        Args:
            n (int): number of dimensions
        """
        reference = _create_tensor(n)
        test_data = xml.load(self.ref_dir + 'tensor{}.xml'.format(n))
        assert np.array_equal(test_data, reference)

    def test_load_arrayofindex(self):
        """Load reference XML file for ARTS type ArrayOfIndex."""
        reference = [1., 2., 3.]
        test_data = xml.load(self.ref_dir + 'arrayofindex.xml')
        assert np.array_equal(test_data, reference)

    def test_load_arrayofindex_binary(self):
        """Load reference binary XML file for ARTS type ArrayOfIndex."""
        reference = [1., 2., 3.]
        test_data = xml.load(self.ref_dir + 'arrayofindex-bin.xml')
        assert np.array_equal(test_data, reference)

    def test_load_arrayofstring(self):
        """Load reference XML file for ARTS type ArrayOfString."""
        reference = ['a', 'bb', 'ccc']
        test_data = xml.load(self.ref_dir + 'arrayofstring.xml')
        assert np.array_equal(test_data, reference)

    def test_load_arrayofvector(self):
        """Load reference XML file for ARTS type ArrayOfVector."""
        reference = [np.arange(1), np.arange(1)]
        test_data = xml.load(self.ref_dir + 'arrayofvector.xml')
        assert np.array_equal(test_data, reference)

    def test_load_comment(self):
        """Load reference XML file storing only a comment."""
        test_data = xml.load(self.ref_dir + 'comment.xml')
        assert test_data is None

    def test_load_arrayofindex_with_comment(self):
        """Load reference XML file for ARTS type ArrayOfIndex with comment."""
        reference = [1., 2., 3.]
        test_data = xml.load(self.ref_dir + 'arrayofindex-comment.xml')
        assert np.array_equal(test_data, reference)


class TestSave:
    """Testing the ARTS XML saving functions.

    This class provides functions to test the saving of XML files. Data is
    created and stored to a temporay file. Afterwards the file is read and
    the data gets compared to the initial data.

    Notes:
        The functions setUp() and tearDown() are run automatically before every
        other function.

    """
    def setup_method(self):
        """Create a temporary file."""
        fd, self.f = mkstemp()
        os.close(fd)

    def teardown_method(self):
        """Delete temporary file."""
        for f in [self.f, self.f + '.gz', self.f + '.bin']:
            if os.path.isfile(f):
                os.remove(f)

    @pytest.mark.parametrize("inttype", (int, np.int32, np.int64))
    def test_save_index(self, inttype):
        """Save Index to file, read it and compare the results."""
        reference = inttype(0)
        xml.save(reference, self.f)
        test_data = xml.load(self.f)
        assert test_data == reference

    def test_save_vector(self):
        """Save Vector to file, read it and compare the results."""
        reference = _create_tensor(1)
        xml.save(reference, self.f)
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    def test_save_vector_binary(self):
        """Save Vector to binary file, read it and compare the results."""
        reference = _create_tensor(1)
        xml.save(reference, self.f, format='binary')
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    def test_save_complex_vector(self):
        """Save complex Vector to file, read it and compare the results."""
        reference = _create_complex_tensor(1)
        xml.save(reference, self.f)
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    def test_save_complex_vector_binary(self):
        """Save complex Vector to file, read it and compare the results."""
        reference = _create_complex_tensor(1)
        xml.save(reference, self.f, format='binary')
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    def test_save_empty_vector(self):
        """Save empty Vector to file, read it and compare the results."""
        reference = _create_empty_tensor(1)
        xml.save(reference, self.f)
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    def test_save_empty_matrix(self):
        """Save empty Matrix to file, read it and compare the results."""
        reference = _create_empty_tensor(2)
        xml.save(reference, self.f)
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    def test_save_matrix(self):
        """Save Matrix to file, read it and compare the results."""
        reference = _create_tensor(2)
        xml.save(reference, self.f)
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    def test_save_matrix_binary(self):
        """Save Matrix to file, read it and compare the results."""
        reference = _create_tensor(2)
        xml.save(reference, self.f, format='binary')
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    @pytest.mark.parametrize('fileformat', ['ascii', 'binary'])
    def test_sparse(self, fileformat):
        """Save Sparse to file, read it and compare the result."""
        reference = _create_sparse(10)
        xml.save(reference, self.f, format=fileformat)
        test_data = xml.load(self.f)
        assert np.array_equal(test_data.toarray(), reference.toarray())

    def test_save_complex_matrix(self):
        """Save complex Matrix to file, read it and compare the results."""
        reference = _create_complex_tensor(2)
        xml.save(reference, self.f)
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    def test_save_complex_matrix_binary(self):
        """Save complex Matrix to file, read it and compare the results."""
        reference = _create_complex_tensor(2)
        xml.save(reference, self.f, format='binary')
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    @pytest.mark.parametrize("inttype", (int, np.int32, np.int64))
    def test_save_arrayofindex(self, inttype):
        """Save ArrayOfIndex to file, read it and compare the results."""
        reference = [inttype(i) for i in [1., 2., 3.]]
        xml.save(reference, self.f)
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    @pytest.mark.parametrize("inttype", (int, np.int32, np.int64))
    def test_save_arrayofindex_binary(self, inttype):
        """Save ArrayOfIndex to binary file, read it and compare the result."""
        reference = [inttype(i) for i in [1., 2., 3.]]
        xml.save(reference, self.f, format='binary')
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    def test_save_arrayofstring(self):
        """Save ArrayOfString to file, read it and compare the results."""
        reference = ['a', 'bb', 'ccc']
        xml.save(reference, self.f)
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    def test_save_arrayofvector(self):
        """Save ArrayOfIndex to file, read it and compare the results."""
        reference = [np.arange(1), np.arange(1)]
        xml.save(reference, self.f)
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    def test_save_gzip(self):
        """Test writing/reading of gzipped files."""
        f = self.f + '.gz'
        ref = np.arange(10)

        xml.save(ref, f)

        assert np.array_equal(ref, xml.load(f))

    def test_save_binary_gzip(self):
        """Check for exception when attempting to write zipped binary file."""
        f = self.f + '.gz'
        ref = np.arange(10)

        with pytest.raises(Exception):
            xml.save(ref, f, format='binary')

    @pytest.mark.parametrize('n', range(3, 8))
    @pytest.mark.parametrize('fileformat', ['ascii', 'binary'])
    def test_save_load_tensor(self, n, fileformat):
        """Save tensor of dimension n to file, read it and compare data to
        reference.

        Args:
            n (int): number of dimensions
            fileformat (str): 'ascii' or 'binary'.

        """
        reference = _create_tensor(n)
        xml.save(reference, self.f, format=fileformat)
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)

    @pytest.mark.parametrize('n', range(3, 8))
    def test_save_empty_tensor(self, n):
        """Save empty tensor of dimension n to file, read it and compare data
        to reference.

        Args:
            n (int): number of dimensions

        """
        reference = _create_empty_tensor(n)
        xml.save(reference, self.f)
        test_data = xml.load(self.f)
        assert np.array_equal(test_data, reference)
