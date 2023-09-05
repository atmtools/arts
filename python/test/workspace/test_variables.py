"""
Test accessing and transferring workspace variables.
"""
import os
import numpy as np
import pytest
import scipy as sp
import pyarts
from pyarts.workspace import Workspace
from pyarts.workspace import global_data as global_data
from pyarts.xml import load, save

class TestVariables:
    """
    Tests the manipulation of workspace variables.
    """
    def setup_method(self):
        """
        This ensures a new Workspace for every test.
        """
        self.dir = os.path.dirname(os.path.realpath(__file__))
        self.ws  = Workspace()
        self.setup_workspace()

    def setup_workspace(self):
        ws = self.ws

        ws.f_grid = 183.0e9 * np.ones(1)

        ws.sensor_los = 180.0 * np.ones((1, 1))
        ws.sensor_pos = 830e3 * np.ones((1, 1))
        ws.sensorOff()

    def test_index_transfer(self):
        """
        Create and set Index WSV.
        """
        self.ws.index_variable = pyarts.arts.Index()
        i = np.random.randint(0, 100)
        self.ws.index_variable = i
        assert self.ws.index_variable == i

    def test_string_transfer(self):
        """
        Create and set String WSV.
        """
        self.ws.string_variable = pyarts.arts.String("string_variable")
        s = "some random string."
        self.ws.string_variable = s
        assert self.ws.string_variable == s

    def test_array_of_index_transfer(self):
        """
        Create and set ArrayOfIndex WSV.
        """
        self.ws.array_of_index_variable = pyarts.arts.ArrayOfIndex()
        i = [np.random.randint(0, 100) for j in range(10)]

        self.ws.array_of_index_variable = i
        assert all(np.array(self.ws.array_of_index_variable.value) == i)

        self.ws.array_of_index_variable = []
        assert all(np.array(self.ws.array_of_index_variable.value) == [])

    def test_array_of_vector_transfer(self):
        """
        Create and set ArrayOfVector WSV.

        """
        self.ws.array_of_vector_variable = pyarts.arts.ArrayOfVector()
        aov = pyarts.xml.load(os.path.join(self.dir, "../xml/reference/arrayofvector.xml"))
        self.ws.array_of_vector_variable = aov
        assert all(np.array(self.ws.array_of_vector_variable) == aov)

    def test_vector_transfer(self):
        """
        Create and set Vector WSV.
        """
        self.ws.vector_variable = pyarts.arts.Vector()
        v = np.random.rand(10)
        self.ws.vector_variable = v
        assert all(self.ws.vector_variable == v)

    def test_matrix_transfer(self):
        """
        Create and set Matrix WSV.
        """
        self.ws.matrix_variable = pyarts.arts.Matrix()
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
                            offsets=[2, 1, 0, -1, -2]).tocsc()
        self.ws.sensor_response = m
        assert np.all(m.toarray() == self.ws.sensor_response.value.toarray())

    def test_time(self):
        """
        Create and set Time variable.
        """
        import datetime
        times = ["2020-01-02 03:04:05", datetime.datetime.now()]
        self.ws.time_1 = pyarts.arts.ArrayOfTime()
        self.ws.ArrayOfTimeNLinSpace(self.ws.time_1, 5, times[0], str(times[1]))

        assert self.ws.time_1[0].time.year == 2020
        assert self.ws.time_1[0].time.month == 1
        assert self.ws.time_1[0].time.day == 2
        assert self.ws.time_1[0].time.hour == 3
        assert self.ws.time_1[0].time.minute == 4
        
        time = pyarts.arts.Time(times[1])
        assert time.time == times[1]

    def test_creation(self):
        """
        Test creation of WSVs.
        """
        self.ws.array_of_index = pyarts.arts.ArrayOfIndex()
        self.ws.array_of_index = pyarts.arts.ArrayOfIndex()
        with pytest.raises(Exception):
            self.ws.array_of_index = pyarts.arts.Numeric()

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

        v = global_data.convert("Index", 1.2)
        assert(v == 1)

        v = global_data.convert("String", "string")
        assert(v == "string")

        v = global_data.convert("Numeric", 1)
        print(type(v))
        assert(type(v) == np.float64)

        v = global_data.convert("Vector", 1.0)
        assert(v.shape == (1,))

        v = global_data.convert("Matrix", 1.0)
        assert(v.shape == (1, 1))

        v = global_data.convert("Tensor3", 1.0)
        assert(v.shape == (1, 1, 1))

        v = global_data.convert("Tensor6", 1.0)
        assert(v.shape == (1, 1, 1, 1, 1, 1))

        v = global_data.convert("ArrayOfArrayOfIndex", 1.0)
        assert(type(v) == list)
        assert(type(v[0]) == list)
        assert(type(v[0][0]) == int)

        v = global_data.convert("ArrayOfArrayOfIndex", 1)
        assert(type(v) == list)
        assert(type(v[0]) == list)
        assert(type(v[0][0]) == int)

    def test_typeerror(self):
        with pytest.raises(TypeError):
            self.ws.f_grid = 1

        with pytest.raises(TypeError):
            self.ws.f_grid = 'c'


if __name__ == "__main__":
    ta = TestVariables()
    ta.setup_method()
    ta.test_creation()
