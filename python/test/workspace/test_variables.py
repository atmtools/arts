"""
Test accessing and transferring workspace variables.
"""
import os
import numpy as np
import pytest
import scipy as sp
import arts
from arts.workspace import Workspace, WorkspaceVariable
from arts.xml import load, save

class TestVariables:
    """
    Tests the manipulation of workspace variables.
    """
    def setup_method(self):
        """
        This ensures a new Workspace for every test.
        """
        self.dir = os.path.dirname(os.path.realpath(__file__))
        self.ws  = Workspace(verbosity = 0)
        self.setup_workspace()

    def setup_workspace(self):
        ws = self.ws
        ws.atmosphere_dim = 1
        ws.p_grid = np.linspace(1e5, 1e3, 21)
        ws.Touch(ws.lat_grid)
        ws.Touch(ws.lon_grid)

        ws.f_grid = 183.0e9 * np.ones(1)
        ws.stokes_dim = 1

        ws.sensor_los = 180.0 * np.ones((1, 1))
        ws.sensor_pos = 830e3 * np.ones((1, 1))
        ws.sensorOff()

    def test_index_transfer(self):
        """
        Create and set Index WSV.
        """
        self.ws.IndexCreate("index_variable")
        i = np.random.randint(0, 100)
        self.ws.index_variable = i
        assert self.ws.index_variable.value == i

    def test_string_transfer(self):
        """
        Create and set String WSV.
        """
        self.ws.StringCreate("string_variable")
        s = "some random string."
        self.ws.string_variable = s
        assert self.ws.string_variable.value == s

    def test_array_of_index_transfer(self):
        """
        Create and set ArrayOfIndex WSV.
        """
        self.ws.ArrayOfIndexCreate("array_of_index_variable")
        i = [np.random.randint(0, 100) for j in range(10)]

        self.ws.array_of_index_variable = i
        assert self.ws.array_of_index_variable.value == i

        self.ws.array_of_index_variable = []
        assert self.ws.array_of_index_variable.value == []

    def test_array_of_vector_transfer(self):
        """
        Create and set ArrayOfVector WSV.

        """
        self.ws.ArrayOfVectorCreate("array_of_vector_variable")
        aov = arts.xml.load(os.path.join(self.dir, "../xml/reference/arrayofvector.xml"))
        self.ws.array_of_vector_variable = aov
        assert self.ws.array_of_vector_variable.value == aov

    def test_vector_transfer(self):
        """
        Create and set Vector WSV.
        """
        self.ws.VectorCreate("vector_variable")
        v = np.random.rand(10)
        self.ws.vector_variable = v
        assert all(self.ws.vector_variable.value == v)

    def test_matrix_transfer(self):
        """
        Create and set Matrix WSV.
        """
        self.ws.MatrixCreate("matrix_variable")
        m = np.random.rand(10, 10)
        self.ws.matrix_variable = m
        assert all(self.ws.matrix_variable.value.ravel() == m.ravel())

    def test_sparse_transfer(self):
        """
        Create and set Sparse WSV.
        """
        n = 100
        d2 = np.ones(n - 2)
        d1 = np.ones(n - 1)
        d = np.ones(n)
        m = sp.sparse.diags(diagonals=[d2, d1, d, d1, d2],
                            offsets=[2, 1, 0, -1, -2])
        self.ws.sensor_response = m
        assert np.all(m.toarray() == self.ws.sensor_response.value.toarray())

    def test_tensor_3(self):
        """
        Create and set Tensor3 variable.
        """
        t_0 = np.random.rand(*([3] * 3))
        self.ws.Tensor3Create("tensor_3")
        self.ws.tensor_3 = t_0
        assert np.all(t_0 == self.ws.tensor_3.value)

    def test_tensor_4(self):
        """
        Create and set Tensor4 variable.
        """
        t_0 = np.random.rand(*([3] * 4))
        t_1 = self.ws.Tensor4Create("tensor_4")
        self.ws.tensor_4 = t_0
        assert np.all(t_0 == self.ws.tensor_4.value)

    def test_tensor_5(self):
        """
        Create and set Tensor5 variable.
        """
        t_0 = np.random.rand(*([3] * 5))
        t_1 = self.ws.Tensor5Create("tensor_5")
        self.ws.tensor_5 = t_0
        assert np.all(t_0 == self.ws.tensor_5.value)

    def test_tensor_6(self):
        """
        Create and set Tensor6 variable.
        """
        t_0 = np.random.rand(*([3] * 6))
        t_1 = self.ws.Tensor6Create("tensor_6")
        self.ws.tensor_6 = t_0
        assert np.all(t_0 == self.ws.tensor_6.value)

    def test_tensor_7(self):
        """
        Create and set Tensor7 variable.
        """
        t_0 = np.random.rand(*([3] * 7))
        self.ws.Tensor7Create("tensor_7")
        self.ws.tensor_7 = t_0
        assert np.all(t_0 == self.ws.tensor_7.value)

    def test_creation(self):
        """
        Test creation of WSVs.
        """
        self.ws.ArrayOfIndexCreate("array_of_index")
        self.ws.ArrayOfIndexCreate("array_of_index")
        with pytest.raises(Exception):
            self.ws.VectorCreate("array_of_index")

    def test_covariance_matrix(self):
        """
        Test manipulation of CorvarianceMatrix objects.
        """
        ws = self.ws

        ws.jacobianInit()
        ws.jacobianAddAbsSpecies(species = "O3",
                                 g1 = ws.p_grid,
                                 g2 = ws.lat_grid,
                                 g3 = ws.lon_grid)
        ws.jacobianAddAbsSpecies(species = "H2O",
                                 g1 = ws.p_grid,
                                 g2 = ws.lat_grid,
                                 g3 = ws.lon_grid)
        ws.jacobianClose()

        ws.covmatDiagonal(out = ws.covmat_block,
                          out_inverse = ws.covmat_block,
                          vars = 10.0 * np.ones(ws.p_grid.value.size))
        ws.covmat_sxAddBlock(block = ws.covmat_block)
        ws.covmatDiagonal(out = ws.covmat_block,
                          out_inverse = ws.covmat_block,
                          vars = 20.0 * np.ones(ws.p_grid.value.size))
        ws.covmat_sxAddBlock(block = ws.covmat_block)

    def test_variable_creation(self):
        """
        Test creation of named and unnambed WSVs.
        """

        # Unnamed variable
        wsv = self.ws.create_variable("Matrix", None)
        self.ws.__setattr__(wsv.name, np.eye(5))
        assert np.all(np.isclose(np.eye(5),
                                 self.ws.__getattr__(wsv.name).value))

        # Named variable
        wsv = self.ws.create_variable("Matrix", "matrix_wsv")
        self.ws.matrix_wsv = np.eye(5)
        assert np.all(np.isclose(np.eye(5), self.ws.matrix_wsv.value))

    def test_variable_set_empty(self):
        """
        Test initialization of workspace variables.
        """
        self.ws.f_grid = np.array([94e9])
        self.ws.f_grid = []
        assert self.ws.f_grid.value.size == 0

    def test_convert(self):
        """
        Test automatic conversion of Python types.
        """

        v = WorkspaceVariable.convert("Index", 1.2)
        assert(v == 1)

        v = WorkspaceVariable.convert("String", "string")
        assert(v == "string")

        v = WorkspaceVariable.convert("Numeric", 1)
        assert(type(v) == np.float64)

        v = WorkspaceVariable.convert("Vector", 1.0)
        assert(v.shape == (1,))

        v = WorkspaceVariable.convert("Matrix", 1.0)
        assert(v.shape == (1, 1))

        v = WorkspaceVariable.convert("Tensor3", 1.0)
        assert(v.shape == (1, 1, 1))

        v = WorkspaceVariable.convert("Tensor6", 1.0)
        assert(v.shape == (1, 1, 1, 1, 1, 1))

        v = WorkspaceVariable.convert("ArrayOfArrayOfIndex", 1.0)
        assert(type(v) == list)
        assert(type(v[0]) == list)
        assert(type(v[0][0]) == int)

        v = WorkspaceVariable.convert("ArrayOfArrayOfIndex", 1)
        return v
