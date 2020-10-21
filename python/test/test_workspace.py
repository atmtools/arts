# -*- encoding: utf-8 -*-
import os

import numpy as np
import pytest
import scipy as sp

import pyarts
from pyarts.workspace import Workspace, arts_agenda
from pyarts.workspace.variables import WorkspaceVariable

def agenda(ws):
    ws.Print(ws.y, 0)

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

    def test_name_collision(self):
        self.ws.VectorSetConstant(self.ws.f_grid, 10, 1.)
        self.ws.VectorCreate("np")
        f_grid = self.ws.f_grid.value
        assert np.all(np.isclose(f_grid, np.ones(10)))
