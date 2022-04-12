"""
Test handling of agendas of the Python interface.
"""
import os
import numpy as np
import pytest
import scipy as sp
import pyarts
from pyarts.workspace import Workspace, arts_agenda
from pyarts.xml import load, save

@arts_agenda
def ppath_agenda(ws):
      ws.Ignore(ws.rte_pos2)
      ws.ppathStepByStep()


class TestAgendas:
    """
    Tests the calling of ARTS workspace methods.
    """
    def setup_method(self):
        """
        This ensures a new Workspace for every test.
        """
        self.ws = Workspace(verbosity = 0)
        self.setup_workspace()

    def setup_workspace(self):
        self.ws.execute_controlfile("artscomponents/clearsky/TestClearSky.arts")

    def test_assignment(self):
        """
        Test assignment of agendas.
        """
        ws = self.ws
        ws.ppath_agenda = ppath_agenda

    def test_include(self):
        ws = self.ws

        @arts_agenda
        def ppath_agenda_inc(ws):
            INCLUDE(ppath_agenda)

        ws.ppath_agenda = ppath_agenda_inc

    def test_execution(self):
        """
        Test definition and execution of agendas.
        """

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

    def test_callback(self):
        """
        Test callbacks by re-implementing iy_space_agenda in Python and
        comparing results of yCalc.
        """
        z_ppath = []

        ws = self.ws

        ws.yCalc()
        y_old = np.copy(ws.y.value)

        import scipy.constants as c

        @arts_agenda(allow_callbacks=True)
        def space_agenda(ws):
            # Since everything happens in Python we need
            # to tell ARTS that we are using all in and outputs.
            ws.Ignore(ws.f_grid)
            ws.Ignore(ws.rtp_pos)
            ws.Ignore(ws.rtp_los)
            ws.Touch(ws.iy)

            # Temperatures and frequency
            t = 2.735
            f = ws.f_grid.value

            # Compute radiances
            c1 = 2.0 * c.h / c.c ** 2
            c2 = c.h / c.k
            b = c1 * f ** 3 / (np.exp(c2 * f / t) - 1.0)

            # Put into iy vector.
            ws.iy = np.zeros((f.size, ws.stokes_dim.value))
            ws.iy.value[:, 0] = b

        # Copy ppath_agenda into workspace.
        ws.iy_space_agenda = space_agenda
        ws.yCalc()

        y_new = np.copy(ws.y.value)

        assert(np.allclose(y_new, y_old))


    def test_callback_2(self):
        """
        Test a very complicated Python callback.
        """

        @arts_agenda(allow_callbacks=True)
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

    def test_unknown_wsv(self):
        """
        Ensure that an exception is thrown when an unknown WSV is
        used inside an agenda.

        This covers https://github.com/atmtools/arts/issues/368
        """
        with pytest.raises(ValueError):
            @arts_agenda
            def my_agenda(ws):
                  ws.UnknownMethod()

    def test_starred(self):
        """
        Test expansion of starred expression.
        """
        @arts_agenda
        def agenda(ws):
            """
            This agenda uses a starred expression.
            """
            ws.IndexSet(*[ws.stokes_dim, 42])

        self.ws.stokes_dim = 0
        agenda.execute(self.ws)
        assert self.ws.stokes_dim.value == 42

    def test_double_starred(self):
        """
        Test expansion of starred expression.
        """
        @arts_agenda
        def agenda(ws):
            """
            This agenda uses a starred expression.
            """
            ws.IndexSet(**{"out" : ws.stokes_dim,
                           "value" : 42})

        self.ws.stokes_dim = 0
        agenda.execute(self.ws)
        assert self.ws.stokes_dim.value == 42

    def test_exception(self):
        """
        Ensure that exception is thrown when a agenda
        variable is set to an invalid value.
        """
        @arts_agenda(allow_callbacks=True)
        def abs_xsec_agenda(ws):
              pass
        self.ws = pyarts.workspace.Workspace()
        with pytest.raises(Exception):
              self.ws.abs_xsec_agenda = abs_xsec_agenda
