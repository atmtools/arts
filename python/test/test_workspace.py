# -*- encoding: utf-8 -*-
import os

import numpy as np
import pytest
import scipy as sp

import arts

try:
    from arts.workspace import Workspace, arts_agenda
    from arts.workspace.variables import WorkspaceVariable
except ImportError:
    skip_arts_tests = True
else:
    skip_arts_tests = False


def agenda(ws):
    ws.Print(ws.y, 0)

@pytest.mark.skipif(skip_arts_tests, reason='ARTS library not available')
class TestWorkspace:
    def setup_method(self):
        """This ensures a new Workspace for every test."""
        self.dir = os.path.dirname(os.path.realpath(__file__))
        self.ws  = Workspace()
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
        self.ws.IndexCreate("index_variable")
        i = np.random.randint(0, 100)
        self.ws.index_variable = i
        assert self.ws.index_variable.value == i

    def test_array_of_index_transfer(self):
        self.ws.ArrayOfIndexCreate("array_of_index_variable")
        i = [np.random.randint(0, 100) for j in range(10)]

        self.ws.array_of_index_variable = i
        assert self.ws.array_of_index_variable.value == i

        self.ws.array_of_index_variable = []
        assert self.ws.array_of_index_variable.value == []

    def test_array_of_vector_transfer(self):
        self.ws.ArrayOfVectorCreate("array_of_vector_variable")
        aov = arts.xml.load(os.path.join(self.dir,
                                                "xml/reference/arrayofvector.xml"))
        self.ws.array_of_vector_variable = aov
        assert self.ws.array_of_vector_variable.value == aov

    def test_string_transfer(self):
        self.ws.StringCreate("string_variable")
        s = "some random string."
        self.ws.string_variable = s
        assert self.ws.string_variable.value == s

    def test_vector_transfer(self):
        self.ws.VectorCreate("vector_variable")
        v = np.random.rand(10)
        self.ws.vector_variable = v
        assert all(self.ws.vector_variable.value == v)

    def test_matrix_transfer(self):
        self.ws.MatrixCreate("matrix_variable")
        m = np.random.rand(10, 10)
        self.ws.matrix_variable = m
        assert all(self.ws.matrix_variable.value.ravel() == m.ravel())

    def test_sparse_transfer(self):
        n = 100
        d2 = np.ones(n - 2)
        d1 = np.ones(n - 1)
        d = np.ones(n)
        m = sp.sparse.diags(diagonals=[d2, d1, d, d1, d2],
                            offsets=[2, 1, 0, -1, -2])
        self.ws.sensor_response = m
        assert np.all(m.todense() == self.ws.sensor_response.value.todense())

    def test_supergeneric_overload_resolution(self):
        self.ws.ArrayOfIndexCreate("array_of_index")
        self.ws.ArrayOfArrayOfIndexCreate("array_of_array_of_index")
        self.ws.array_of_index = [1, 2, 3]
        self.ws.Append(self.ws.array_of_array_of_index, self.ws.array_of_index)
        self.ws.Append(self.ws.array_of_array_of_index, self.ws.array_of_index)

    def test_creation(self):
        self.ws.ArrayOfIndexCreate("array_of_index")
        self.ws.ArrayOfIndexCreate("array_of_index")
        with pytest.raises(Exception):
            self.ws.VectorCreate("array_of_index")

    def test_wsm_error(self):
        with pytest.raises(Exception):
            self.ws.yCalc()

    def test_doc(self):
        repr(self.ws.yCalc)

    def test_agenda(self):

        self.ws.atmosphere_dim = 1

        @arts_agenda
        def add_1(ws):
            ws.IndexAdd(ws.atmosphere_dim,
                        ws.atmosphere_dim,
                        1)
        add_1.execute(self.ws)

        assert self.ws.atmosphere_dim.value == 2

        add_1.append(add_1)
        add_1.execute(self.ws)

        assert self.ws.atmosphere_dim.value == 4

        args = [self.ws.atmosphere_dim, self.ws.atmosphere_dim, 1]

        @arts_agenda
        def add_2(ws):
            ws.IndexAdd(*args)

        add_2.execute(self.ws)

        assert self.ws.atmosphere_dim.value == 5


    def test_execute_controlfile(self):

        dir = os.path.dirname(os.path.realpath(__file__))
        test_dir = os.path.join(dir, "test_files")
        self.ws.WriteXML("ascii", np.array([1.0]),
                         os.path.join(test_dir, "vector.xml"))
        os.chdir(test_dir)
        self.ws.execute_controlfile("controlfile.arts")

        os.remove(os.path.join(test_dir, "vector.xml"))

    def test_supergeneric_overload_failure(self):
        with pytest.raises(Exception):
            self.ws.NumericCreate("numeric_wsv")
            self.ws.StringCreate("string_wsv")
            self.ws.Copy(self.ws.string_wsv, self.ws.numeric_wsv)

    def test_tensor_3(self):
        t_0 = np.random.rand(*([3] * 3))
        self.ws.Tensor3Create("tensor_3")
        self.ws.tensor_3 = t_0
        assert np.all(t_0 == self.ws.tensor_3.value)

    def test_tensor_4(self):
        t_0 = np.random.rand(*([3] * 4))
        t_1 = self.ws.Tensor4Create("tensor_4")
        self.ws.tensor_4 = t_0
        assert np.all(t_0 == self.ws.tensor_4.value)

    def test_tensor_5(self):
        t_0 = np.random.rand(*([3] * 5))
        t_1 = self.ws.Tensor5Create("tensor_5")
        self.ws.tensor_5 = t_0
        assert np.all(t_0 == self.ws.tensor_5.value)

    def test_tensor_6(self):
        t_0 = np.random.rand(*([3] * 6))
        t_1 = self.ws.Tensor6Create("tensor_6")
        self.ws.tensor_6 = t_0
        assert np.all(t_0 == self.ws.tensor_6.value)

    def test_tensor_7(self):
        t_0 = np.random.rand(*([3] * 7))
        self.ws.Tensor7Create("tensor_7")
        self.ws.tensor_7 = t_0
        assert np.all(t_0 == self.ws.tensor_7.value)

    def test_execute_controlfile(self):

        dir = os.path.dirname(os.path.realpath(__file__))
        test_dir = os.path.join(dir, "test_files")
        self.ws.WriteXML("ascii", np.array([1.0]),
                         os.path.join(test_dir, "vector.xml"))
        os.chdir(test_dir)

        agenda = self.ws.execute_controlfile("controlfile.arts")
        self.ws.foo = "not bar"

        @arts_agenda
        def execute(ws):
            ws.FlagOff(ws.jacobian_do)
            ws.StringSet(ws.foo, "still not bar")
            INCLUDE("controlfile.arts")
            INCLUDE(agenda)

        self.ws.execute_agenda(execute)

        assert self.ws.foo.value == "bar"
        os.remove(os.path.join(test_dir, "vector.xml"))

    def test_covariance_matrix(self):
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



    def test_variable_set_empty(self):
        self.ws.f_grid = np.array([94e9])
        self.ws.f_grid = []
        assert self.ws.f_grid.value.size == 0

    def test_variable_creation(self):

        # Unnamed variable
        wsv = self.ws.create_variable("Matrix", None)
        self.ws.__setattr__(wsv.name, np.eye(5))
        assert np.all(np.isclose(np.eye(5),
                                 self.ws.__getattr__(wsv.name).value))

        # Named variable
        wsv = self.ws.create_variable("Matrix", "matrix_wsv")
        self.ws.matrix_wsv = np.eye(5)
        assert np.all(np.isclose(np.eye(5), self.ws.matrix_wsv.value))

    def test_wsv_setattr(self):
        wsv = self.ws.atmosphere_dim
        wsv.value = 12
        assert self.ws.atmosphere_dim.value == 12

    def test_convert(self):

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

    def test_callbacks(self):

        @arts_agenda
        def agenda(ws):
            """
            This agenda sets a workspace variable in a very
            obscure way.
            """

            class Foo:
                def __init__(self, ws):
                    self.ws = ws

                def ooo(self):
                    self.ws.IndexSet(ws.stokes_dim, 42)

            foo = Foo(ws)
            ws.IndexSet(ws.stokes_dim, 21)
            foo.ooo()

        agenda.execute(self.ws)

        assert self.ws.stokes_dim.value == 42

    def test_contiguous_arrays(self):
        x = np.linspace(0, 1, 256)

        xf = np.asarray(x, order='F')
        self.ws.f_grid = xf
        assert np.array_equal(self.ws.f_grid.value, xf)

        self.ws.f_grid = x[::2]
        assert np.array_equal(self.ws.f_grid.value, x[::2])

        self.ws.f_grid = np.ascontiguousarray(x[::2])
        assert np.array_equal(self.ws.f_grid.value, x[::2])
