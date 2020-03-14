# DEFINITIONS:  -*-sh-*-
#
# Agenda predefinitions for iy_surface_sub_agendaX and surface_rtprop_sub_agendaX
#
# No need to include this file for standard calculations, but handy to include
# if you operate with different surface types. This file then provides dummy
# definitions of all agendas of concern, and you just need to set up the
# agendas you actually use.
#
# The iy_surface_sub_agenda-s are all set to use the corresponding
# surface_rtprop_agenda.
#
# The surface_rtprop_sub_agenda-s all contain just Ignore and Touch calls, and
# will trigger an error if used.
#
# Authors: Patrick Eriksson

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)


@arts_agenda
def iy_surface_sub_agenda0(ws):
    ws.AgendaExecute(ws.surface_rtprop_sub_agenda0)
    ws.iySurfaceRtpropCalc()


ws.iy_surface_sub_agenda0 = iy_surface_sub_agenda0


@arts_agenda
def iy_surface_sub_agenda1(ws):
    ws.AgendaExecute(ws.surface_rtprop_sub_agenda1)
    ws.iySurfaceRtpropCalc()


ws.iy_surface_sub_agenda1 = iy_surface_sub_agenda1


@arts_agenda
def iy_surface_sub_agenda2(ws):
    ws.AgendaExecute(ws.surface_rtprop_sub_agenda2)
    ws.iySurfaceRtpropCalc()


ws.iy_surface_sub_agenda2 = iy_surface_sub_agenda2


@arts_agenda
def iy_surface_sub_agenda3(ws):
    ws.AgendaExecute(ws.surface_rtprop_sub_agenda3)
    ws.iySurfaceRtpropCalc()


ws.iy_surface_sub_agenda3 = iy_surface_sub_agenda3


@arts_agenda
def iy_surface_sub_agenda4(ws):
    ws.AgendaExecute(ws.surface_rtprop_sub_agenda4)
    ws.iySurfaceRtpropCalc()


ws.iy_surface_sub_agenda4 = iy_surface_sub_agenda4


@arts_agenda
def iy_surface_sub_agenda5(ws):
    ws.AgendaExecute(ws.surface_rtprop_sub_agenda5)
    ws.iySurfaceRtpropCalc()


ws.iy_surface_sub_agenda5 = iy_surface_sub_agenda5


@arts_agenda
def surface_rtprop_sub_agenda0(ws):
    ws.Error("This agenda just contain dummy calls, you need to define it properly.")
    ws.Ignore(ws.f_grid)
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)
    ws.Ignore(ws.surface_type_aux)
    ws.Touch(ws.surface_emission)
    ws.Touch(ws.surface_los)
    ws.Touch(ws.surface_rmatrix)


ws.surface_rtprop_sub_agenda0 = surface_rtprop_sub_agenda0

ws.Copy(ws.surface_rtprop_sub_agenda1, ws.surface_rtprop_sub_agenda0)
ws.Copy(ws.surface_rtprop_sub_agenda2, ws.surface_rtprop_sub_agenda0)
ws.Copy(ws.surface_rtprop_sub_agenda3, ws.surface_rtprop_sub_agenda0)
ws.Copy(ws.surface_rtprop_sub_agenda4, ws.surface_rtprop_sub_agenda0)
ws.Copy(ws.surface_rtprop_sub_agenda5, ws.surface_rtprop_sub_agenda0)
