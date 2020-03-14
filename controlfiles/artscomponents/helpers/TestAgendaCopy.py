#
# Testing and demonstrating how to use user-defined agenda variants.
#
# Jana Mendrok 2013-02-26

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.AgendaCreate("propmat_clearsky_agenda__OnTheFly")


@arts_agenda
def propmat_clearsky_agenda__OnTheFly(ws):
    ws.Ignore(ws.rtp_mag)
    ws.Ignore(ws.rtp_los)
    ws.Ignore(ws.rtp_nlte)
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddOnTheFly()


ws.propmat_clearsky_agenda__OnTheFly = propmat_clearsky_agenda__OnTheFly

ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
ws.AgendaCreate("mybatch")


@arts_agenda
def mybatch(ws):
    ws.Touch(ws.y)
    ws.Touch(ws.y_aux)
    ws.Touch(ws.jacobian)
    ws.Print(ws.ybatch_index)


ws.mybatch = mybatch

ws.IndexSet(ws.ybatch_start, 0)
ws.IndexSet(ws.ybatch_n, 10)
ws.Copy(ws.ybatch_calc_agenda, ws.mybatch)
ws.ybatchCalc()
