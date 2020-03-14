#
# Testing whether (and demonstrating how) nested forloops work.
# Result: Amazing, but YES, it works! :-)
#
# Jana Mendrok 2013-02-26

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.StringCreate("out")
ws.IndexCreate("outlevel")
ws.IndexSet(ws.outlevel, 0)
ws.AgendaCreate("forloop_agenda_inner")


@arts_agenda
def forloop_agenda_inner(ws):
    ws.StringSet(ws.out, "inner")
    ws.Print(ws.out, ws.outlevel)
    ws.Print(ws.forloop_index, ws.outlevel)


ws.forloop_agenda_inner = forloop_agenda_inner

ws.AgendaCreate("forloop_agenda_outer")


@arts_agenda
def forloop_agenda_outer(ws):
    ws.StringSet(ws.out, "outer")
    ws.Print(ws.out, ws.outlevel)
    ws.Print(ws.forloop_index, ws.outlevel)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_inner)
    ws.ForLoop(ws.forloop_agenda, 0, 9, 3)


ws.forloop_agenda_outer = forloop_agenda_outer

ws.Copy(ws.forloop_agenda, ws.forloop_agenda_outer)
ws.ForLoop(ws.forloop_agenda, 0, 2, 1)
