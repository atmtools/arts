# DEFINITIONS:  -*-sh-*-
#
############
# Agenda predefinitions for surface_rtprop_agenda
#
# This is a sub-file to agendas.arts, separated for better overview.
#
############
#
# Authors: Jana Mendrok
#

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
##################################
#    surface_rtprop_agenda
##################################
#####
# blackbody surface
#####
# Fastest option, but only a good choice only if troposphere is opaque at the
# frequency in question.
# - surface skin temperature interpolated from t_surface
#
ws.AgendaCreate("surface_rtprop_agenda__Blackbody_SurfTFromt_surface")


@arts_agenda
def surface_rtprop_agenda__Blackbody_SurfTFromt_surface(ws):
    ws.InterpSurfaceFieldToPosition(out=ws.surface_skin_t, field=ws.t_surface)
    ws.surfaceBlackbody()


ws.surface_rtprop_agenda__Blackbody_SurfTFromt_surface = (
    surface_rtprop_agenda__Blackbody_SurfTFromt_surface
)

# - surface skin temperature interpolated from atmospheric t_field
#
ws.AgendaCreate("surface_rtprop_agenda__Blackbody_SurfTFromt_field")


@arts_agenda
def surface_rtprop_agenda__Blackbody_SurfTFromt_field(ws):
    ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
    ws.surfaceBlackbody()


ws.surface_rtprop_agenda__Blackbody_SurfTFromt_field = (
    surface_rtprop_agenda__Blackbody_SurfTFromt_field
)

#####
# specular reflecting surfaces (no polarization)
#####
# - preset surface reflectivity
# - surface skin temperature interpolated from t_surface
#
ws.AgendaCreate("surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface")


@arts_agenda
def surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface(ws):
    ws.specular_losCalc()
    ws.InterpSurfaceFieldToPosition(out=ws.surface_skin_t, field=ws.t_surface)
    ws.surfaceFlatScalarReflectivity()


ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface = (
    surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface
)

# - preset surface reflectivity
# - surface skin temperature interpolated from atmospheric t_field
#
ws.AgendaCreate("surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_field")


@arts_agenda
def surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_field(ws):
    ws.specular_losCalc()
    ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
    ws.surfaceFlatScalarReflectivity()


ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_field = (
    surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_field
)

#####
# specular reflecting surfaces (considering polarization)
#####
# - preset surface reflectivity
# - surface skin temperature interpolated from t_surface
#
ws.AgendaCreate("surface_rtprop_agenda__Specular_WithPol_ReflFix_SurfTFromt_surface")


@arts_agenda
def surface_rtprop_agenda__Specular_WithPol_ReflFix_SurfTFromt_surface(ws):
    ws.specular_losCalc()
    ws.InterpSurfaceFieldToPosition(out=ws.surface_skin_t, field=ws.t_surface)
    ws.surfaceFlatReflectivity()


ws.surface_rtprop_agenda__Specular_WithPol_ReflFix_SurfTFromt_surface = (
    surface_rtprop_agenda__Specular_WithPol_ReflFix_SurfTFromt_surface
)

#####
# Lambertian reflecting surfaces
#####
# - preset surface reflectivity
# - surface skin temperature interpolated from t_surface
# - surface normal derived via specular_losCalc (when set manually, specular_losCalc can be removed here)
# - default number of angles, the scattering function will be evaluated on
#
ws.AgendaCreate("surface_rtprop_agenda__lambertian_ReflFix_SurfTFromt_surface")


@arts_agenda
def surface_rtprop_agenda__lambertian_ReflFix_SurfTFromt_surface(ws):
    ws.InterpSurfaceFieldToPosition(out=ws.surface_skin_t, field=ws.t_surface)
    ws.specular_losCalc()
    ws.surfaceLambertianSimple()


ws.surface_rtprop_agenda__lambertian_ReflFix_SurfTFromt_surface = (
    surface_rtprop_agenda__lambertian_ReflFix_SurfTFromt_surface
)

# - preset surface reflectivity
# - surface skin temperature interpolated from atmospheric t_field
# - surface normal derived via specular_losCalc (when set manually, specular_losCalc can be removed here)
# - default number of angles, the scattering function will be evaluated on
#
ws.AgendaCreate("surface_rtprop_agenda__lambertian_ReflFix_SurfTFromt_field")


@arts_agenda
def surface_rtprop_agenda__lambertian_ReflFix_SurfTFromt_field(ws):
    ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
    ws.specular_losCalc()
    ws.surfaceLambertianSimple()


ws.surface_rtprop_agenda__lambertian_ReflFix_SurfTFromt_field = (
    surface_rtprop_agenda__lambertian_ReflFix_SurfTFromt_field
)
