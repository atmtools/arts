# DEFINITIONS:  -*-sh-*-
#
# Test and compares three manners to set spectral_radiance_field
#
# 2019-10-27, Patrick Eriksson

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
# on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
ws.NumericSet(ws.ppath_lmax, 100.0)
# Geometrical path calculation (i.e., refraction neglected)
#
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# Standard RT agendas
#
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
# Definition of species
#
ws.abs_speciesSet(species=["N2-SelfContStandardType", "O2-PWR98", "H2O-PWR98"])
# No line data needed here
#
ws.abs_lines_per_speciesSetEmpty()
# Atmosphere
#
ws.AtmosphereSet1D()
ws.VectorNLogSpace(ws.p_grid, 161, 101300.0, 1.0)
ws.AtmRawRead(basename="testdata/tropical")
#
ws.AtmFieldsCalc()
# Define (a blacbody) surface
#
ws.MatrixSet(ws.z_surface, array([[0.0]]))
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Blackbody_SurfTFromt_field)
# Map this to variables used by DISORT
#
ws.VectorSet(ws.surface_scalar_reflectivity, array([0.0]))
ws.VectorExtractFromMatrix(ws.rtp_pos, ws.z_surface, 0, "row")
ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
# Frequencies and Stokes dim.
#
ws.IndexSet(ws.stokes_dim, 1)
ws.VectorSet(ws.f_grid, array([1.84e11]))
# Stuff not used
#
ws.jacobianOff()
ws.sensorOff()
# Create a za_grid
#
ws.AngularGridsSetFluxCalc(N_za_grid=20)
# Perform checks
#
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
ws.lbl_checkedCalc()
ws.FlagOn(ws.sensor_checked)
# We start with pure clear-sky calculations
#
ws.cloudboxOff()
ws.cloudbox_checkedCalc()
#
ws.Tensor3Create("dummy")
#
ws.spectral_radiance_fieldClearskyPlaneParallel(trans_field=ws.dummy)
#
ws.Tensor7Create("ref_field")
ws.Copy(ws.ref_field, ws.spectral_radiance_field)
# WriteXML( "binary", z_field, "z.xml" )
# WriteXML( "binary", spectral_radiance_field, "f1.xml" )
# Set spectral_radiance_field as full DISORT calculation
#
ws.DisortCalcClearsky(nstreams=8)
# WriteXML( "binary", spectral_radiance_field, "f2.xml" )
ws.Compare(ws.ref_field, ws.spectral_radiance_field, 2e-18)
# Copy( ref_field, spectral_radiance_field )
# For last test, we create an empty cloudbox up to about 10 km
#
ws.cloudboxSetManually(
    p1=1100000.0, p2=22000.0, lat1=-1.0, lat2=-1.0, lon1=-1.0, lon2=-1.0
)
ws.pnd_fieldZero()
ws.scat_data_checkedCalc()
# Run DISORT and use a third method
#
ws.DisortCalc(nstreams=8)
ws.Copy(ws.iy_cloudbox_agenda, ws.iy_cloudbox_agenda__QuadInterpField1D)
ws.spectral_radiance_fieldExpandCloudboxField()
# WriteXML( "binary", spectral_radiance_field, "f3.xml" )
ws.Compare(ws.ref_field, ws.spectral_radiance_field, 2e-18)
