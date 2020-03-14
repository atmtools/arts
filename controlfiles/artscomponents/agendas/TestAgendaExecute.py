# Make sure scoping of variables works correctly for
# AgendaExecute and AgendaExecuteExclusive

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)


@arts_agenda
def test_agenda(ws):
    ws.AgendaExecute(ws.g0_agenda)


ws.test_agenda = test_agenda


@arts_agenda
def g0_agenda(ws):
    ws.NumericSet(ws.lat, -1.0)
    ws.NumericSet(ws.g0, 20.0)
    ws.Ignore(ws.lat)
    ws.Ignore(ws.lon)


ws.g0_agenda = g0_agenda

ws.NumericSet(ws.lat, 0.0)
ws.Touch(ws.lon)
ws.NumericCreate("ref")
ws.NumericSet(ws.g0, 10.0)
ws.AgendaExecute(ws.test_agenda)
ws.NumericSet(ws.ref, 10.0)
ws.Compare(ws.g0, ws.ref, 0.0, "Error in Agenda Scoping")
ws.AgendaExecute(ws.g0_agenda)
ws.NumericSet(ws.ref, 20.0)
ws.Compare(ws.g0, ws.ref, 0.0, "Error in Agenda Scoping")
ws.NumericSet(ws.g0, 10.0)
ws.AgendaExecuteExclusive(ws.test_agenda)
ws.NumericSet(ws.ref, 10.0)
ws.Compare(ws.g0, ws.ref, 0.0, "Error in Agenda Scoping")
ws.AgendaExecuteExclusive(ws.g0_agenda)
ws.NumericSet(ws.ref, 20.0)
ws.Compare(ws.g0, ws.ref, 0.0, "Error in Agenda Scoping")
ws.AgendaExecuteExclusive(ws.g0_agenda)
ws.NumericSet(ws.ref, 0.0)
ws.Compare(ws.lat, ws.ref, 0.0, "Error in Agenda Scoping")
ws.MatrixCreate("mref")
ws.MatrixSet(ws.mref, np.array([[4.0, 5.0, 6.0]]))
ws.VectorSet(ws.x, [])
ws.IndexSet(ws.jacobian_do, 0)
ws.IndexSet(ws.inversion_iteration_counter, 0)
ws.MatrixSet(ws.jacobian, np.array([[1.0, 2.0, 3.0]]))


@arts_agenda
def inversion_iterate_agenda(ws):
    ws.Ignore(ws.x)
    ws.Ignore(ws.jacobian_do)
    ws.Ignore(ws.jacobian)
    ws.Ignore(ws.inversion_iteration_counter)
    ws.VectorSet(ws.yf, np.array([4.0, 5.0, 6.0]))
    ws.MatrixSet(ws.jacobian, np.array([[4.0, 5.0, 6.0]]))


ws.inversion_iterate_agenda = inversion_iterate_agenda

ws.AgendaExecute(ws.inversion_iterate_agenda)
ws.Compare(ws.jacobian, ws.mref, 0.0, "Error in Agenda InOut Scoping")
ws.Print(ws.jacobian, 0)
