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
# This case is for Mars and specifically tests
#
#  - surface altitude:
#      - reading in raw 3D data, preprocess into z_surface variable within a 3D
#         global case (!but no atmfields_checkedCalc - as, on a global
#         scale this fails since altitude reference of altitude data and
#         atmospheric data is inconsistent. with plenty of negative altitudes in
#         surface altitude data and atm data starting only at 5m ARTS rejects
#         the data due to gap between atmo and surface!) (CASE A-1)
#      - reading in raw 3D data, preprocess into z_surface variable within a 3D
#         regional case containing only surface altitudes, and perform
#         atmfields_checkedCalc (CASE A-2, B)
#  - surface refractive index data (CASE B, C):
#      - reading 3D raw data
#      - in surface_rt_prop_agenda deriving surface emission/reflection field
#         from complex refractive index data (using surfaceFlatRefractiveIndex)
#      - all in 3D only for regional (CASE B) and global (CASEs C) cases
#  - surface temperature (CASE C):
#      - reading in raw 3D data, preprocess into t_surface variable, and within
#        surface_rt_prop_agenda derive surface_skin_t from t_surface
#      - loop over all scenarios (4seasons x 2daytimes x 3dustloads)
#      - all in global 3D
#
# Jana Mendrok 2013-02-26

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_mars.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# the 3D geo grid for global case
ws.VectorCreate("lat_grid3D_glob")
ws.VectorCreate("lon_grid3D_glob")
ws.VectorLinSpace(ws.lat_grid3D_glob, -90.0, 90.0, 5.0)
ws.VectorLinSpace(ws.lon_grid3D_glob, 0.0, 360.0, 10.0)
# the 3D geo grid for regional case (z_surface>5m)
ws.VectorCreate("lat_grid3D_reg")
ws.VectorCreate("lon_grid3D_reg")
ws.VectorLinSpace(ws.lat_grid3D_reg, -49.0, -18.0, 1.0)
ws.VectorLinSpace(ws.lon_grid3D_reg, 208.0, 298.0, 2.0)
ws.GriddedField2Create("gf2tmp")
ws.MatrixCreate("mtmp")
ws.IndexCreate("itmp")
ws.IndexCreate("ncases")
ws.StringCreate("basename")
ws.StringSet(ws.basename, "planets/Mars/MPS/")
ws.StringCreate("zsurfname")
ws.StringSet(ws.zsurfname, "Mars.z_surface")
ws.StringCreate("risurfname")
ws.StringSet(ws.risurfname, "Mars.surface_complex_refr_index_field")
ws.StringCreate("atmcase")
ws.StringCreate("caseext")
ws.StringCreate("casefull")
ws.StringCreate("tsurfname")
ws.StringSet(ws.tsurfname, ".t_surface")
ws.StringCreate("atmext")
ws.StringSet(ws.atmext, ".sol-avg")
# Arrays with (sub)case names
ws.ArrayOfStringCreate("seasoncasearray")
ws.ArrayOfStringSet(
    ws.seasoncasearray, ["Mars.Ls0", "Mars.Ls90", "Mars.Ls180", "Mars.Ls270"]
)
ws.ArrayOfStringCreate("timecasearray")
ws.ArrayOfStringSet(ws.timecasearray, [".day", ".night"])
ws.ArrayOfStringCreate("dustcasearray")
ws.ArrayOfStringSet(ws.dustcasearray, [".dust-high", ".dust-low", ".dust-medium"])
# a vector for holding reference RT results
ws.VectorCreate("yREFERENCE")
# and some further settings in order to be able to do an RT calc
#####
ws.jacobianOff()
ws.cloudboxOff()
ws.IndexSet(ws.stokes_dim, 1)
ws.VectorSet(ws.f_grid, array([3.0e11]))
ws.sensorOff()
ws.StringSet(ws.iy_unit, "PlanckBT")
# we manually select a minumim set of basic atm data (main atm constituents)
ws.abs_speciesSet(species=["CO2"])
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
ws.MatrixSetConstant(ws.mtmp, ws.itmp, 1, -19.0)
ws.Append(ws.sensor_pos, ws.mtmp, "trailing")
# sensor longitutde
ws.MatrixSetConstant(ws.mtmp, ws.itmp, 1, 210.0)
ws.Append(ws.sensor_pos, ws.mtmp, "trailing")
#####
# CASE A
# applying surface altitude. fixed surface reflectivity, t_surface from t_field.
# first part global, then reducing to region with z_surface>5m
#####
# some stuff to get a basic atmosphere
#####
ws.Extract(ws.casefull, ws.seasoncasearray, 0)
ws.Copy(ws.atmcase, ws.casefull)
ws.Extract(ws.casefull, ws.timecasearray, 0)
ws.Append(ws.atmcase, ws.casefull)
ws.Extract(ws.casefull, ws.dustcasearray, 0)
ws.Append(ws.atmcase, ws.casefull)
ws.Copy(ws.casefull, ws.basename)
ws.Append(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, "/")
ws.Append(ws.casefull, ws.caseext)
ws.Append(ws.atmcase, ws.atmext)
ws.Append(ws.casefull, ws.atmcase)
ws.Append(ws.casefull, ws.caseext)
ws.Append(ws.casefull, ws.atmcase)
# Print( casefull, 0 )
ws.AtmRawRead(basename=ws.casefull)
ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
ws.AtmosphereSet3D()
#####
# CASE A-1
#####
ws.Copy(ws.lat_grid, ws.lat_grid3D_glob)
ws.Copy(ws.lon_grid, ws.lon_grid3D_glob)
ws.AtmFieldsCalcExpand1D()
# reading the surface altitude field
ws.Copy(ws.casefull, ws.basename)
ws.Append(ws.casefull, ws.zsurfname)
ws.ReadXML(ws.gf2tmp, ws.casefull)
ws.GriddedFieldLatLonRegrid(ws.gf2tmp, ws.lat_grid, ws.lon_grid, ws.gf2tmp)
ws.FieldFromGriddedField(ws.z_surface, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf2tmp)
#####
# CASE A-2
#####
ws.Copy(ws.lat_grid, ws.lat_grid3D_reg)
ws.Copy(ws.lon_grid, ws.lon_grid3D_reg)
ws.AtmFieldsCalcExpand1D()
# reading the surface altitude field
ws.Copy(ws.casefull, ws.basename)
ws.Append(ws.casefull, ws.zsurfname)
ws.ReadXML(ws.gf2tmp, ws.casefull)
ws.GriddedFieldLatLonRegrid(ws.gf2tmp, ws.lat_grid, ws.lon_grid, ws.gf2tmp)
ws.FieldFromGriddedField(ws.z_surface, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf2tmp)
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
# now we need to do some RT calc in order to APPLY the reflectivity data
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.VectorSet(ws.surface_scalar_reflectivity, array([0.4]))
# surface temp from atmospheric t_field
@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
    # Print( surface_skin_t, 0 )
    ws.surfaceFlatScalarReflectivity()


ws.surface_rtprop_agenda = surface_rtprop_agenda

ws.yCalc()
# Print( y, 0 )
#####
# CASE B
# applying surface altitude and refractive index. t_surface from t_field.
# still regional only.
#####
# reading surface refractive index data (GriddedField5). no regridding here. we
#  need to do that in form of 1D only inside the surface_rtprop_agenda, as we
#  only know then which exact point(s) we need. for better traceability (and
#  since here this isn't just a temporary field), we create a dedicated
#  workspace variable for this data.
ws.Copy(ws.casefull, ws.basename)
ws.Append(ws.casefull, ws.risurfname)
ws.GriddedField5Create("ri_surface")
ws.ReadXML(ws.ri_surface, ws.casefull)


@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
    # Print( surface_skin_t, 0 )
    ws.Select(ws.lat_true, ws.rtp_pos, [1])
    ws.Select(ws.lon_true, ws.rtp_pos, [2])
    ws.surface_complex_refr_indexFromGriddedField5(
        complex_refr_index_field=ws.ri_surface
    )
    ws.surfaceFlatRefractiveIndex()


ws.surface_rtprop_agenda = surface_rtprop_agenda

ws.yCalc()
# Print( y, 0 )
# WriteXML( in=y )
ws.ReadXML(out=ws.yREFERENCE, filename="TestSurf_Mars.y.xml")
ws.Compare(ws.y, ws.yREFERENCE, 1e-06)
#####
# CASE C
# applying refractive index and t_surface. z_surface from lowest atmo-z.
# returning to global case
#####
ws.AtmosphereSet3D()
ws.Copy(ws.lat_grid, ws.lat_grid3D_glob)
ws.Copy(ws.lon_grid, ws.lon_grid3D_glob)


@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.InterpSurfaceFieldToPosition(out=ws.surface_skin_t, field=ws.t_surface)
    # Print( surface_skin_t, 0 )
    ws.Select(ws.lat_true, ws.rtp_pos, [1])
    ws.Select(ws.lon_true, ws.rtp_pos, [2])
    ws.surface_complex_refr_indexFromGriddedField5(
        complex_refr_index_field=ws.ri_surface
    )
    ws.surfaceFlatRefractiveIndex()


ws.surface_rtprop_agenda = surface_rtprop_agenda

# we go with several nested foorloop through the different cases.
ws.AgendaCreate("forloop_agenda_dust")


@arts_agenda
def forloop_agenda_dust(ws):
    # construct atmcase name III (Mars.LsXX.YY.dust-ZZ)
    ws.Extract(ws.casefull, ws.dustcasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    # keep the casestring till dust and make upper-level folder name
    ws.Append(ws.basename, ws.atmcase)
    ws.StringSet(ws.caseext, "/")
    ws.Append(ws.basename, ws.caseext)
    ws.Copy(ws.casefull, ws.basename)
    ws.Append(ws.casefull, ws.atmcase)
    # reading the surface temperature field
    ws.Append(ws.casefull, ws.tsurfname)
    ws.ReadXML(ws.gf2tmp, ws.casefull)
    ws.GriddedFieldLatLonRegrid(ws.gf2tmp, ws.lat_grid, ws.lon_grid, ws.gf2tmp)
    # reading atmospheric field data
    ws.Append(ws.atmcase, ws.atmext)
    ws.Append(ws.basename, ws.atmcase)
    ws.Append(ws.basename, ws.caseext)
    ws.Append(ws.basename, ws.atmcase)
    ws.AtmRawRead(basename=ws.basename)
    ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
    ws.AtmFieldsCalcExpand1D()
    ws.Extract(ws.z_surface, ws.z_field, 0)
    ws.FieldFromGriddedField(
        ws.t_surface, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf2tmp
    )
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.sensor_checkedCalc()
    # Print( casefull, 0 )
    ws.yCalc()
    # Print( y, 0 )


ws.forloop_agenda_dust = forloop_agenda_dust

ws.AgendaCreate("forloop_agenda_time")


@arts_agenda
def forloop_agenda_time(ws):
    # construct atmcase name II (Mars.LsXX.d/n)
    ws.Extract(ws.casefull, ws.timecasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_dust)
    ws.nelemGet(ws.ncases, ws.dustcasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_time = forloop_agenda_time

ws.AgendaCreate("forloop_agenda_season")


@arts_agenda
def forloop_agenda_season(ws):
    # construct atmcase name I (Mars.LsXX)
    ws.Extract(ws.casefull, ws.seasoncasearray, ws.forloop_index)
    ws.Copy(ws.atmcase, ws.casefull)
    ws.Copy(ws.forloop_agenda, ws.forloop_agenda_time)
    ws.nelemGet(ws.ncases, ws.timecasearray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.forloop_agenda_season = forloop_agenda_season

ws.nelemGet(ws.ncases, ws.seasoncasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.forloop_agenda_season)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
