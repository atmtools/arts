# --------------------------------------------------------------------
# Test batch calculations for AVHRR. See also TestHIRS.arts.
# Gerrit Holl, July 2011
# Based on TestHIRS.arts
# --------------------------------------------------------------------

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
# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# standard surface agenda (i.e., make use of surface_rtprop_agenda)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# blackbody surface with skin temperature interpolated from t_surface field
ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Blackbody_SurfTFromt_field)
ws.StringCreate("satellite")
ws.ArrayOfIndexCreate("channels")
ws.ArrayOfIndexCreate("views")
ws.StringCreate("hitran_file")
ws.NumericCreate("f_grid_spacing")
# Select here which satellite you want
# ---
ws.StringSet(ws.satellite, "NOAA19")
# Select here which channels you want
# ---
#
# ARTS 0 --> AVHRR 3B
# ARTS 1 --> AVHRR 4
# ARTS 2 --> AVHRR 5
ws.ArrayOfIndexSet(ws.channels, [1, 2])
# Select here which views you want.
# AVHRR LAC/FRAC has 1024 views
# ---
ws.ArrayOfIndexSet(ws.views, [0, 256, 512, 768, 1023])
ws.StringSet(ws.hitran_file, "/storage3/data/catalogue/hitran/hitran2004/HITRAN04.par")
# Set frequency grid spacing
# (See comments in hirs_reference.arts concerning useful values here)
# ---
ws.NumericSet(ws.f_grid_spacing, 500000000.0)
# Basic settings (already needed in sensor part)
# ---
# This example assumes 1D
ws.AtmosphereSet1D()
# scalar RT
ws.IndexSet(ws.stokes_dim, 1)
ws.execute_controlfile("avhrr_reference.arts")
# Set up absorption
# =================
ws.ReadXML(ws.batch_atm_fields_compact, "testdata/garand_profiles.xml.gz")
# add constant profiles for O2 and N2
ws.batch_atm_fields_compactAddConstant(name="abs_species-O2", value=0.2095)
ws.batch_atm_fields_compactAddConstant(name="abs_species-N2", value=0.7808)
# Setup lookup table
# ---
ws.abs_lookupSetupBatch()
ws.abs_xsec_agenda_checkedCalc()
ws.abs_lookupCalc()
# Set propmat_clearsky_agenda to use lookup table
# ---
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# Set up RT calculation
# =====================
# Set the agenda for batch calculations:
# ---
#
@arts_agenda
def ybatch_calc_agenda(ws):
    # Extract the atmospheric profiles for this case:
    ws.Extract(ws.atm_fields_compact, ws.batch_atm_fields_compact, ws.ybatch_index)
    ws.AtmFieldsAndParticleBulkPropFieldFromCompact()
    # get some surface properties from corresponding atmospheric fields
    ws.Extract(ws.z_surface, ws.z_field, 0)
    ws.Extract(ws.t_surface, ws.t_field, 0)
    # No scattering
    ws.cloudboxOff()
    # No jacobian calculations
    ws.jacobianOff()
    # Perform RT calculations
    # ---
    ws.atmfields_checkedCalc()
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.sensor_checkedCalc()
    ws.yCalc()
    # Convert the measurement from radiance units to Planck Tb:
    ws.StringSet(ws.iy_unit, "PlanckBT")
    ws.yApplyUnit()


ws.ybatch_calc_agenda = ybatch_calc_agenda

# Set number of batch cases (uncomment according to your needs):
# this one for all cases in atmosphere batch
# nelemGet( ybatch_n, batch_atm_fields_compact )
# for testing, ALL cases take too long. so we do only the first 2 (still takes
# ~15min on an 8-node system.
ws.IndexSet(ws.ybatch_n, 2)
# Execute the batch calculations:
# ---
ws.propmat_clearsky_agenda_checkedCalc()
ws.ybatchCalc()
# Store result matrix:
# ---
ws.WriteXML("ascii", ws.ybatch)
