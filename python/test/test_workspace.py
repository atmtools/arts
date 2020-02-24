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
