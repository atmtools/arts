# DEFINITIONS:  -*-sh-*-
#
# ARTS controlfile testing consistency of use of built-in and read-from-file
# isotopologue ratios. It calculates radiances with applying 3 different
# approaches of getting & using isotopologue ratios (IR):
# (a) the classical case: using ARTS built-in IR
# (b) using identical IR as from (a), but read in from file
# (c) for H2O-162 (=HDO), we use IR_new=IR_orig*2 and VMR_new=VMR_orig. As in
#     ARTS absorption = VMR*IR*LineStrength*LineShapeFct, this should give
#     (largely) identical results to (a) and (b).
#
# 2012-09-14 Jana Mendrok

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
ws.VectorLinSpace(ws.f_grid, 488000000000.0, 491000000000.0, 10000000.0)
ws.VectorNLogSpace(ws.p_grid, 41, 100000.0, 1.0)
# monchromatic pencilbeam 1D scalar clearsky RT without Jacobians
ws.AtmosphereSet1D()
ws.IndexSet(ws.stokes_dim, 1)
ws.jacobianOff()
ws.cloudboxOff()
ws.sensorOff()
# output radiances in Planck-Tb
ws.StringSet(ws.iy_unit, "RJBT")
# sensor position and LOS
ws.IndexCreate("n_tan")
ws.IndexSet(ws.n_tan, 5)
ws.VectorCreate("z_tan")
ws.VectorNLinSpace(ws.z_tan, ws.n_tan, 50000.0, 10000.0)
# Create vector of zenith angles for selected tangent altitudes
ws.MatrixSetConstant(ws.sensor_pos, ws.n_tan, 1, 600000.0)
ws.VectorCreate("za")
ws.VectorZtanToZa1D(ws.za, ws.sensor_pos, ws.refellipsoid, ws.atmosphere_dim, ws.z_tan)
ws.Matrix1ColFromVector(ws.sensor_los, ws.za)
# On-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# Blackbody surface
ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Blackbody_SurfTFromt_field)
# Definition of species: here we neglect continua
ws.abs_speciesSet(species=["H2O-161", "H2O-162", "N2", "O2"])
# use planetary catalogue
ws.ReadSplitARTSCAT(basename="spectroscopy/Perrin/", fmin=0.0, fmax=3000000000000.0)
ws.abs_linesSetCutoff(option="ByLine", value=750000000000.0)
ws.abs_linesSetNormalization(option="VVH")
ws.abs_lines_per_speciesCreateFromLines()
#######
# CASE (a): H2O only with builtIn setup and standard profiles
#######
# Atmospheric scenario
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalc()
# Surface at bottom of defined atosphere
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.lbl_checkedCalc()
# check basic atm&surface setup of
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
# now do the RT calc
ws.isotopologue_ratiosInitFromBuiltin()
# WriteXML("ascii", isotopologue_ratios, "IRfromARTSBuiltin.xml")
ws.yCalc()
# WriteXML( "ascii", y, "y.IsoRatios_Builtin.xml" )
ws.VectorCreate("yA")
ws.Copy(ws.yA, ws.y)
#######
# CASE (b): H2O with (unchanged) IR from file (relevant H2O only) and AtmFiles explicitly read
#######
# we need some container for vmr reading
ws.GriddedField3Create("vmr1_raw")
ws.ArrayOfGriddedField3Create("vmrall_raw")
# Read Atmospheric scenario files
ws.ReadXML(ws.t_field_raw, "testdata/tropical.t.xml")
ws.ReadXML(ws.z_field_raw, "testdata/tropical.z.xml")
# twice read&append H2O as we use default profile for both defined H2O tag groups
ws.ReadXML(ws.vmr1_raw, "testdata/tropical.H2O.xml")
ws.Append(ws.vmrall_raw, ws.vmr1_raw)
ws.ReadXML(ws.vmr1_raw, "testdata/tropical.H2O.xml")
ws.Append(ws.vmrall_raw, ws.vmr1_raw)
ws.ReadXML(ws.vmr1_raw, "testdata/tropical.N2.xml")
ws.Append(ws.vmrall_raw, ws.vmr1_raw)
ws.ReadXML(ws.vmr1_raw, "testdata/tropical.O2.xml")
ws.Append(ws.vmrall_raw, ws.vmr1_raw)
ws.Copy(ws.vmr_field_raw, ws.vmrall_raw)
# calc fields from raw gridded fields
ws.AtmFieldsCalc()
# check basic atm&surface setup of
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
# now do the RT calc
ws.ReadXML(ws.isotopologue_ratios, "IsoRatios_H2Oorig.xml")
# WriteXML("ascii", isotopologue_ratios, "IRfromH2Oorig.xml")
ws.yCalc()
# WriteXML( "ascii", y, "y.IsoRatios_H2Oorig.xml" )
ws.VectorCreate("yB")
ws.Copy(ws.yB, ws.y)
ws.Compare(ws.yA, ws.yB, 1e-20)
#######
# CASE (c): H2O-162 with modified IR (*2) from file and VMR profile (/2)
#######
ws.ArrayOfGriddedField3Create("vmrall_raw2")
# Read Atmospheric scenario files
ws.ReadXML(ws.t_field_raw, "testdata/tropical.t.xml")
ws.ReadXML(ws.z_field_raw, "testdata/tropical.z.xml")
# here read separate files for H2O and HDO tag groups
ws.ReadXML(ws.vmr1_raw, "testdata/tropical.H2O.xml")
ws.Append(ws.vmrall_raw2, ws.vmr1_raw)
ws.ReadXML(ws.vmr1_raw, "tropical.H2O_50%.xml")
ws.Append(ws.vmrall_raw2, ws.vmr1_raw)
ws.ReadXML(ws.vmr1_raw, "testdata/tropical.N2.xml")
ws.Append(ws.vmrall_raw2, ws.vmr1_raw)
ws.ReadXML(ws.vmr1_raw, "testdata/tropical.O2.xml")
ws.Append(ws.vmrall_raw2, ws.vmr1_raw)
ws.Copy(ws.vmr_field_raw, ws.vmrall_raw2)
# calc fields from raw gridded fields
ws.AtmFieldsCalc()
# check basic atm&surface setup of
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
# now do the RT calc
ws.ReadXML(ws.isotopologue_ratios, "IsoRatios_H2Omani.xml")
ws.yCalc()
# WriteXML( "ascii", y, "y.IsoRatios_H2Omani.xml" )
ws.VectorCreate("yC")
ws.Copy(ws.yC, ws.y)
ws.Compare(ws.yA, ws.yC, 0.01)
