# -*- coding: utf-8 -*-
"""Testing the covariance matrix class.
"""
import shutil

import pytest
import numpy as np
import scipy as sp
import os
from tempfile import mkstemp
from pyarts.covariancematrix import Block, CovarianceMatrix
from pyarts.xml import load, save

class TestCovarianceMatrix:

    def setup_method(self):
        # Temporary file
        fd, self.f = mkstemp()
        os.close(fd)

        # Simple covariance matrix for testing
        b1 = Block(0, 0, 0, 0, False, np.random.normal(size = (10, 10)))
        b2 = Block(1, 1, 10, 10, False, sp.sparse.identity(10))
        self.covmat = CovarianceMatrix([b1, b2])

    def test_xml_io(self):
        save(self.covmat, self.f)
        covmat2 = load(self.f)

        def compare_matrices(args):
            b1, b2 = args
            m1 = b1.matrix
            m2 = b2.matrix
            if isinstance(m1, sp.sparse.spmatrix):
                m1 = m1.todense()
                m2 = m2.todense()
                print(m1)
            return np.allclose(m1, m2)

        assert(all(map(compare_matrices, zip(self.covmat.blocks, covmat2.blocks))))

    def test_to_dense(self):
        m = self.covmat.to_dense()
        assert(np.allclose(m[:10, :10], self.covmat.blocks[0].matrix))
        assert(np.allclose(m[10:, 10:], self.covmat.blocks[1].matrix.toarray()))

    def teardown_method(self):
        # Remove temp file
        os.remove(self.f)
