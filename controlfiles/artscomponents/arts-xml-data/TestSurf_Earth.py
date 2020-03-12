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
# This case is for Earth and specifically tests
#
#  - surface altitude: reading in raw 3D data, preprocess into z_surface variable
#      - within a 3D (global) case (CASE A)
#      - as 1D (extract one point) case (CASE B)
#  - land-sea-mask: used as basic data to derive land/sea dependent
#     reflectivities (maps mask value to refl=0.05 for land and 0.4 for sea)
#      - reading 3D raw data
#      - applying surface_rt_prop agenda on-the run extraction of mask value and
#         derivation of reflectivity at surface reflection point (CASE A-2;
#         CASE A-1 using fixed surface reflectivity for comparison)
#      - extracting at single point & used in 1D case (CASE B; sensor viewing
#        and lat/lon_true adjusted such that similar region observed as for
#        CASEs A)
#
#  - CASEs A-1, A-2, and B RT calculation results compared for consistency
#     (setup such that A-1 and A-2 have the exactly same sensor setup, A-2
#     looking to ocean and A-1 reflectivity set to the one of ocean; B has same
#     viewing geometry looking at similar point (but one point only!) as
#     cases A. Hence, results of A-1 and A-2 shall be identical, while results
#     from B shall be "close".)
#
# Jana Mendrok 2013-02-26

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# parameters for converting land/sea mask to surface reflectivities
#  mask has land=1, sea=0
# later on we will convert those values to refl_land and refl_sea by
#  refl=a*maskval+b, where a=refl_land-refl_sea and b=refl_sea)
ws.NumericCreate("refl_land")
ws.NumericCreate("refl_sea")
ws.NumericSet(ws.refl_land, 0.05)
# emiss=0.95
ws.NumericSet(ws.refl_sea, 0.4)
# emiss=0.60
# the 3D geo grid to test
ws.VectorCreate("lat_grid3D")
ws.VectorCreate("lon_grid3D")
ws.VectorLinSpace(ws.lat_grid3D, -90.0, 90.0, 18.0)
ws.VectorLinSpace(ws.lon_grid3D, 0.0, 360.0, 18.0)
# VectorLinSpace( lon_grid3D, 20, 90, 2 )
# alternatively, the geocoordinates for the 1D case (before switching on, check
#  to not overwrite them further down)
# VectorSet( lat_true, [32.] )
# VectorSet( lon_true, [40.] )
ws.GriddedField2Create("gf2tmp")
ws.MatrixCreate("mtmp")
ws.VectorCreate("vtmp")
ws.NumericCreate("ntmp")
ws.IndexCreate("itmp")
ws.StringCreate("atmcase")
ws.StringSet(ws.atmcase, "planets/Earth/Fascod/tropical/tropical")
# StringSet( atmcase, "planets/Earth/Fascod/midlatitude-summer/midlatitude-summer" )
ws.StringCreate("surfpath")
ws.StringSet(ws.surfpath, "planets/Earth/ECMWF/ERA40/")
ws.StringCreate("zsurfname")
ws.StringSet(ws.zsurfname, "SurfaceAltitude_ERA40_1.0Degree")
ws.StringCreate("surfmaskname")
ws.StringSet(ws.surfmaskname, "LandSeaMask_ERA40_1.0Degree")
# StringSet( surfmaskname, "LandSeaMask_ERA40_0.25Degree" )
ws.StringCreate("casefull")
ws.StringCreate("caseext")
# some stuff to get a basic atmosphere
#####
# we manually select a minumim set of basic atm data (main atm constituents)
ws.abs_speciesSet(species=["O2", "N2"])
ws.AtmRawRead(basename=ws.atmcase)
ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
# and some further settings in order to be able to do an RT calc
#####
ws.jacobianOff()
ws.cloudboxOff()
ws.IndexSet(ws.stokes_dim, 1)
ws.VectorSet(ws.f_grid, np.array([3.0e11]))
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
# LOS zenith angle
ws.MatrixSet(ws.sensor_los, np.array([[180.0], [130.0], [115.0], [113.8]]))
# LOS azimuth angle
# MatrixSet( mtmp,       [])
ws.nrowsGet(ws.itmp, ws.sensor_los)
ws.MatrixSetConstant(ws.mtmp, ws.itmp, 1, 90.0)
ws.Append(ws.sensor_los, ws.mtmp, "trailing")
# sensor altitude
ws.MatrixSetConstant(ws.sensor_pos, ws.itmp, 1, 600000.0)
# sensor latitude
ws.MatrixSetConstant(ws.mtmp, ws.itmp, 1, 0.0)
ws.Append(ws.sensor_pos, ws.mtmp, "trailing")
# sensor longitutde
ws.MatrixSetConstant(ws.mtmp, ws.itmp, 1, 170.0)
ws.Append(ws.sensor_pos, ws.mtmp, "trailing")
#####
# CASEs A: we start out with 3D (both fields are 3D)
#####
ws.AtmosphereSet3D()
ws.Copy(ws.lat_grid, ws.lat_grid3D)
ws.Copy(ws.lon_grid, ws.lon_grid3D)
ws.AtmFieldsCalcExpand1D(vmr_zeropadding=1)
# reading the surface altitude field
ws.Copy(ws.casefull, ws.surfpath)
ws.Append(ws.casefull, ws.zsurfname)
ws.ReadXML(ws.gf2tmp, ws.casefull)
ws.GriddedFieldLatLonRegrid(ws.gf2tmp, ws.lat_grid, ws.lon_grid, ws.gf2tmp)
ws.FieldFromGriddedField(ws.z_surface, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf2tmp)
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
# reading surface mask; no regridding here. we need to do that as 1D inside the
#  surface_rtprop_agenda, as we only know then which exact point(s) we need. for
#  better traceability (and since here this isn't just a temporary field), we
#  create a dedicated workspace variable for the mask data.
ws.Copy(ws.casefull, ws.surfpath)
ws.Append(ws.casefull, ws.surfmaskname)
ws.GriddedField2Create("lsmask")
ws.ReadXML(ws.lsmask, ws.casefull)
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
# now we need to do some RT calc in order to APPLY the reflectivity data
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
# CASE A-1
#####
ws.VectorSet(ws.surface_scalar_reflectivity, np.array([0.4]))
ws.Copy(
    ws.surface_rtprop_agenda,
    ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_field,
)
ws.yCalc()
# WriteXML( "ascii", y, "TestSurf_Earth.y.3D.fixRef.xml" )
ws.VectorCreate("y_3D_fixRef")
ws.Copy(ws.y_3D_fixRef, ws.y)
# Print( y )
# CASE A-2
#####
@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
    ws.Select(ws.lat_true, ws.rtp_pos, [1])
    ws.Select(ws.lon_true, ws.rtp_pos, [2])
    ws.GriddedFieldLatLonRegrid(ws.lsmask, ws.lat_true, ws.lon_true, ws.lsmask)
    ws.FieldFromGriddedField(ws.mtmp, ws.p_grid, ws.lat_true, ws.lon_true, ws.lsmask)
    ws.VectorExtractFromMatrix(ws.surface_scalar_reflectivity, ws.mtmp, 0, "row")
    ws.Copy(ws.ntmp, ws.refl_sea)
    ws.NumericScale(ws.ntmp, ws.ntmp, -1.0)
    ws.NumericAdd(ws.ntmp, ws.ntmp, ws.refl_land)
    ws.VectorScale(
        ws.surface_scalar_reflectivity, ws.surface_scalar_reflectivity, ws.ntmp
    )
    ws.VectorAddScalar(
        ws.surface_scalar_reflectivity, ws.surface_scalar_reflectivity, ws.refl_sea
    )
    # Print( rte_pos, 0 )
    # Print( surface_scalar_reflectivity, 0 )
    ws.surfaceFlatScalarReflectivity()


ws.surface_rtprop_agenda = surface_rtprop_agenda

ws.yCalc()
# WriteXML( "ascii", y, "TestSurf_Earth.y.3D.extRef.xml" )
ws.VectorCreate("y_3D_extRef")
ws.Copy(ws.y_3D_extRef, ws.y)
# Print( y )
ws.Compare(ws.y_3D_fixRef, ws.y_3D_extRef, 1e-06)
#####
# CASE B: now for the fun of it, we do a 1D case, too
#####
# basic atmospheric fields
ws.AtmosphereSet1D()
ws.AtmFieldsCalc(vmr_zeropadding=1)
# sensor in 1D (zenith angles like for 3D. no azimtuh of course)
ws.VectorExtractFromMatrix(ws.lat_true, ws.sensor_pos, 1, "column")
ws.Select(ws.lat_true, ws.lat_true, [0])
ws.VectorExtractFromMatrix(ws.lon_true, ws.sensor_pos, 2, "column")
ws.Select(ws.lon_true, ws.lon_true, [0])
ws.VectorExtractFromMatrix(ws.vtmp, ws.sensor_los, 0, "column")
ws.Matrix1ColFromVector(ws.sensor_los, ws.vtmp)
ws.nrowsGet(ws.itmp, ws.sensor_los)
ws.MatrixSetConstant(ws.sensor_pos, ws.itmp, ws.atmosphere_dim, 600000.0)
# reading the surface altitude field (need to do that again as we didn't keep
#  the "raw" gridded field) and extracting a single point for 1D
ws.Copy(ws.casefull, ws.surfpath)
ws.Append(ws.casefull, ws.zsurfname)
ws.ReadXML(ws.gf2tmp, ws.casefull)
ws.GriddedFieldLatLonRegrid(ws.gf2tmp, ws.lat_true, ws.lon_true, ws.gf2tmp)
ws.FieldFromGriddedField(ws.z_surface, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf2tmp)
# Print( z_surface )
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
# reading surface mask (need to do that again as we didn't keep the "raw"
#  gridded field) and extracting a single point for 1D
ws.Copy(ws.casefull, ws.surfpath)
ws.Append(ws.casefull, ws.surfmaskname)
ws.ReadXML(ws.gf2tmp, ws.casefull)
ws.GriddedFieldLatLonRegrid(ws.gf2tmp, ws.lat_true, ws.lon_true, ws.gf2tmp)
ws.FieldFromGriddedField(ws.mtmp, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.gf2tmp)
ws.VectorExtractFromMatrix(ws.surface_scalar_reflectivity, ws.mtmp, 0, "row")
# converting the land/sea mask value to surface reflectivity
ws.Copy(ws.ntmp, ws.refl_sea)
ws.NumericScale(ws.ntmp, ws.ntmp, -1.0)
ws.NumericAdd(ws.ntmp, ws.ntmp, ws.refl_land)
ws.VectorScale(ws.surface_scalar_reflectivity, ws.surface_scalar_reflectivity, ws.ntmp)
ws.VectorAddScalar(
    ws.surface_scalar_reflectivity, ws.surface_scalar_reflectivity, ws.refl_sea
)
# Print( surface_scalar_reflectivity )
# now we need to do some RT calc in order to APPLY the reflectivity data
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
# Copy( surface_rtprop_agenda, surface_rtprop_agenda__Blackbody_SurfTFromt_field )
# Copy( surface_rtprop_agenda, surface_rtprop_agenda__lambertian_ReflFix_SurfTFromt_field )
ws.Copy(
    ws.surface_rtprop_agenda,
    ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_field,
)
ws.yCalc()
# WriteXML( "ascii", y, "TestSurf_Earth.y.1D.xml" )
# Print( y )
ws.Compare(ws.y, ws.y_3D_extRef, 0.1)
