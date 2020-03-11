#
# Testing functionality (meeting format requirements, etc.) of surface related
#  data.
#
# General test setup: reading in raw data (including a basic atmosphere),
#  extracting/regridding, executing standard pre-RT calc internal test method
#  atmfields_checkedCalc, and performing some RT simulations in order to apply data
#  that has no dedicated check method (e.g., when only used through agendas like
#  surface reflectivity).
#
#
# This case is for Venus and specifically tests
#
#  - surface temperature (CASE A): reading in raw 3D data, preprocess into
#     t_surface variable, and within surface_rt_prop_agenda derive surface_skin_t
#     from t_surface (CASE A-0 first does without t_surface, taking
#     surface_skin_t from atmospheric t_field). all in 3D.
#  - surface refractive index data (CASE B):
#      - reading 3D raw data
#      - in surface_rt_prop_agenda deriving surface emission/reflection field
#         from complex refractive index data (using surfaceFlatRefractiveIndex)
#      - all in 3D only
#
# Jana Mendrok 2013-02-26

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_venus.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# the 3D geo grid to test
ws.VectorCreate("lat_grid3D")
ws.VectorCreate("lon_grid3D")
ws.VectorLinSpace(ws.lat_grid3D, -90.0, 90.0, 18.0)
ws.VectorLinSpace(ws.lon_grid3D, -20.0, 350.0, 18.0)
ws.GriddedField2Create("gf2tmp")
ws.MatrixCreate("mtmp")
ws.IndexCreate("itmp")
ws.StringCreate("basename")
ws.StringSet(ws.basename, "planets/Venus/MPS/")
ws.StringCreate("atmcase")
ws.StringSet(ws.atmcase, "Venus.spicav.night/Venus.spicav.night")
ws.StringCreate("tsurfname")
ws.StringSet(ws.tsurfname, "Venus.t_surface")
ws.StringCreate("risurfname")
ws.StringSet(ws.risurfname, "Venus.surface_complex_refr_index_field")
ws.StringCreate("casefull")
# some stuff to get a basic atmosphere
#####
# we manually select a minumim set of basic atm data (main atm constituents)
ws.abs_speciesSet(species=["CO2"])
ws.Copy(ws.casefull, ws.basename)
ws.Append(ws.casefull, ws.atmcase)
ws.AtmRawRead(basename=ws.casefull)
ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
# and some further settings in order to be able to do an RT calc
#####
ws.jacobianOff()
ws.cloudboxOff()
ws.IndexSet(ws.stokes_dim, 1)
ws.VectorSet(ws.f_grid, array([3.0e11]))
ws.sensorOff()
ws.StringSet(ws.iy_unit, "PlanckBT")
ws.ReadSplitARTSCAT(basename="spectroscopy/Perrin/", fmin=0.0, fmax=1000000000000.0)
ws.abs_lines_per_speciesCreateFromLines()
# and agenda settings needed for RT calc
#####
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# sensor placed over Maxwell Montes region scanning from the high to low surface
#  RI region
# LOS zenith angle
ws.MatrixSet(ws.sensor_los, array([[180.0], [130.0], [115.0], [113.8]]))
# LOS azimuth angle
# MatrixSet( mtmp,       [])
ws.nrowsGet(ws.itmp, ws.sensor_los)
ws.MatrixSetConstant(ws.mtmp, ws.itmp, 1, 120.0)
ws.Append(ws.sensor_los, ws.mtmp, "trailing")
# sensor altitude
ws.MatrixSetConstant(ws.sensor_pos, ws.itmp, 1, 600000.0)
# sensor latitude
ws.MatrixSetConstant(ws.mtmp, ws.itmp, 1, 60.0)
ws.Append(ws.sensor_pos, ws.mtmp, "trailing")
# sensor longitutde
ws.MatrixSetConstant(ws.mtmp, ws.itmp, 1, 60.0)
ws.Append(ws.sensor_pos, ws.mtmp, "trailing")
#####
# CASEs A: first only using t_surface data. reflectivity fixed.
#####
ws.AtmosphereSet3D()
ws.Copy(ws.lat_grid, ws.lat_grid3D)
ws.Copy(ws.lon_grid, ws.lon_grid3D)
ws.AtmFieldsCalcExpand1D(vmr_zeropadding=1)
ws.Extract(ws.z_surface, ws.z_field, 1)
# reading and regridding the surface temperature field
ws.Copy(ws.casefull, ws.basename)
ws.Append(ws.casefull, ws.tsurfname)
ws.ReadXML(ws.gf2tmp, ws.casefull)
ws.GriddedFieldLatLonRegrid(ws.gf2tmp, ws.lat_grid, ws.lon_grid, ws.gf2tmp)
ws.FieldFromGriddedField(ws.t_surface, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf2tmp)
# reading surface refractive index data (GriddedField4). no regridding here. we
#  need to do that in form of 1D only inside the surface_rtprop_agenda, as we
#  only know then which exact point(s) we need. for better traceability (and
#  since here this isn't just a temporary field), we create a dedicated
#  workspace variable for this data.
ws.Copy(ws.casefull, ws.basename)
ws.Append(ws.casefull, ws.risurfname)
ws.GriddedField5Create("ri_surface")
ws.ReadXML(ws.ri_surface, ws.casefull)
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
# now we need to do some RT calc in order to APPLY the reflectivity data
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
# CASE A-0 (we actually don't use t_surface, but surface_skin_t from t_field
#####
ws.VectorSet(ws.surface_scalar_reflectivity, array([0.4]))
# surface temp from atmospheric t_field
@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
    ws.surfaceFlatScalarReflectivity()


ws.surface_rtprop_agenda = surface_rtprop_agenda

# yCalc
# Print( y )
# CASE A
#####
# surface temp from t_surface field
@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.InterpSurfaceFieldToPosition(out=ws.surface_skin_t, field=ws.t_surface)
    ws.surfaceFlatScalarReflectivity()


ws.surface_rtprop_agenda = surface_rtprop_agenda

ws.yCalc()
# Print( y )
#####
# CASE B: now also applying surface reflection properties from surface complex
#  refractive index data
#####
@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.InterpSurfaceFieldToPosition(out=ws.surface_skin_t, field=ws.t_surface)
    ws.Select(ws.lat_true, ws.rtp_pos, [1])
    ws.Select(ws.lon_true, ws.rtp_pos, [2])
    ws.surface_complex_refr_indexFromGriddedField5(
        complex_refr_index_field=ws.ri_surface
    )
    ws.surfaceFlatRefractiveIndex()


ws.surface_rtprop_agenda = surface_rtprop_agenda

ws.yCalc()
# Print( y )
