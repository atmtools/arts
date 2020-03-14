#
# Check functionality of isotopologue ratio from file implementation and
#  functionality of given data in file
#
# We take an atmospheric scenario from Earth
#   CASE A: do a small RT calculation with this based on ARTS built-in
#           (valid for Earth) isotopologue ratios, write them to file
#   CASE B: re-read the isotopologue ratios written in CASE A and repeat the
#           previous calculation (this should give identical results)
#   CASE C(1-3): read the isotopologue ratios prepared for Mars, Venus, and Jupiter and repeat #           the above calculation again (no result check, only being functional)
#
# Jana Mendrok 2013-02-26

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

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
ws.VectorLinSpace(ws.f_grid, 10000000000.0, 3000000000000.0, 10000000000.0)
# WriteXML( "ascii", f_grid )
ws.VectorNLogSpace(ws.p_grid, 41, 100000.0, 0.1)
# monchromatic pencilbeam 1D scalar clearsky RT without Jacobians
ws.AtmosphereSet1D()
ws.IndexSet(ws.stokes_dim, 1)
ws.jacobianOff()
ws.cloudboxOff()
ws.sensorOff()
# output radiances in Planck-Tb
ws.StringSet(ws.iy_unit, "PlanckBT")
# On-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# Blackbody surface
ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Blackbody_SurfTFromt_field)
# set atmospheric scenario
ws.StringCreate("atmcase")
ws.StringSet(ws.atmcase, "planets/Earth/Fascod/tropical/tropical")
# derive abs species from scenario data
ws.abs_speciesDefineAllInScenario(basename=ws.atmcase)
ws.WriteXML(ws.output_file_format, ws.abs_species, "", 0)
# Atmospheric scenario
ws.AtmRawRead(basename=ws.atmcase)
ws.AtmFieldsCalc()
# Surface at bottom of defined atosphere
ws.Extract(ws.z_surface, ws.z_field, 0)
# use planetary catalogue
ws.ReadSplitARTSCAT(basename="spectroscopy/Perrin/", fmin=0.0, fmax=3000000000000.0)
ws.abs_linesSetCutoff(option="ByLine", value=750000000000.0)
ws.abs_linesSetNormalization(option="VVH")
ws.abs_lines_per_speciesCreateFromLines()
# sensor position and LOS
# set the sensor exactly at the surface (as we are in 1D this has exactly 1 element)
ws.Copy(ws.sensor_pos, ws.z_surface)
# also one LOS direction. slant uplooking.
ws.MatrixSetConstant(ws.sensor_los, 1, 1, 60.0)
# check basic atm&surface setup of
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.lbl_checkedCalc()
#######
# CASE A: BuiltIn (Earth) isotoplogue ratios
#######
# now initialize the isotopologue ratios from builtin
ws.isotopologue_ratiosInitFromBuiltin()
# do the RT calc
# yCalc
ws.Touch(ws.y)
# WriteXML( "ascii", y, "y.IsoRatios_EarthBuiltin.xml" )
ws.VectorCreate("yA")
ws.Copy(ws.yA, ws.y)
# store the isotopologue ratios from builtin such that we can read them in in the next step
ws.WriteXML("ascii", ws.isotopologue_ratios)
#######
# CASE B: isotopologue ratios (Earth) read from file
#######
# read in the isotopologue ratios from file
ws.ReadXML(ws.isotopologue_ratios)
# repeat the RT calc
# yCalc
# WriteXML( "ascii", y, "y.IsoRatios_EarthReadFromFile.xml" )
ws.VectorCreate("yB")
ws.Copy(ws.yB, ws.y)
ws.Compare(ws.yA, ws.yB, 1e-20)
#######
# CASEs C: repeat reading in and doing RT with isotopologue ratios from file for
#  the different planets (isotopologue data from toolbox data collection)
#######
ws.VectorCreate("yREFERENCE")
ws.ReadXML(ws.isotopologue_ratios, "planets/Venus/isotopratio_Venus")
ws.yCalc()
# WriteXML( in=y, filename="TestPlanetIsoRatios.Venus.y_reference.xml" )
ws.ReadXML(out=ws.yREFERENCE, filename="TestPlanetIsoRatios.Venus.y_reference.xml")
ws.Compare(ws.yREFERENCE, ws.y, 2e-09)
ws.ReadXML(ws.isotopologue_ratios, "planets/Mars/isotopratio_Mars")
ws.yCalc()
# WriteXML( in=y, filename="TestPlanetIsoRatios.Mars.y_reference.xml" )
ws.ReadXML(out=ws.yREFERENCE, filename="TestPlanetIsoRatios.Mars.y_reference.xml")
ws.Compare(ws.y, ws.yREFERENCE, 2e-09)
ws.ReadXML(ws.isotopologue_ratios, "planets/Jupiter/isotopratio_Jupiter")
ws.yCalc()
# WriteXML( in=y, filename="TestPlanetIsoRatios.Jupiter.y_reference.xml" )
ws.ReadXML(out=ws.yREFERENCE, filename="TestPlanetIsoRatios.Jupiter.y_reference.xml")
ws.Compare(ws.y, ws.yREFERENCE, 2e-09)
