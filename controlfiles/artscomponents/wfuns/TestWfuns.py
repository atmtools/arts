# DEFINITIONS:  -*-sh-*-
#
# This is a test of weighting function calculations.
#
# The test is based on the Odin-SMR 501 GHz case found in another folder.
#
# Author: Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
##############################################################################
#
# Initial part
#
##############################################################################
# Select frequency band here:
#
# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 1)
# 1D atmosphere
#
ws.AtmosphereSet1D()
ws.jacobianOff()
ws.execute_controlfile("instruments/odinsmr/odinsmr_501.arts")
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
# absorption from LUT
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# Copy( propmat_clearsky_agenda, propmat_clearsky_agenda__OnTheFly )
# ---- Atmospheric scenario --------------------------------------------------
# A pressure grid rougly matching 0 to 80 km in 250 m steps.
#
ws.VectorNLogSpace(ws.p_grid, 321, 100000.0, 1.0)
# Atmospheric profiles here taken from Fascod
ws.AtmRawRead(basename="testdata/tropical")
#
ws.AtmFieldsCalc()
# Get ground altitude (z_surface) from z_field
ws.Extract(ws.z_surface, ws.z_field, 0)
# ---- Create absorption table -----------------------------------------------
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_lines_per_speciesSetNormalization(option="VVH")
ws.abs_lines_per_speciesSetCutoff(option="ByLine", value=750000000000.0)
ws.AbsInputFromAtmFields()
ws.abs_speciesSet(abs_species=ws.abs_nls, species=[])
ws.VectorSet(ws.abs_nls_pert, [])
ws.VectorSet(ws.abs_t_pert, [])
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.abs_lookupCalc()
# ---- Sensor position and LOS -----------------------------------------------
# Number of tangent altitudes
ws.IndexCreate("n_tan")
ws.IndexSet(ws.n_tan, 2)
# Sensor position, with platform altitude set to 600 km
ws.MatrixSetConstant(ws.sensor_pos, ws.n_tan, 1, 600000.0)
# LOS, specified by the corresponding geometrical tangent altitudes
# Tangent altitudes will be equally spaced between 50 and 20 km
ws.VectorCreate("z_tan")
ws.VectorNLinSpace(ws.z_tan, ws.n_tan, 50000.0, 20000.0)
ws.VectorCreate("za")
ws.VectorZtanToZa1D(ws.za, ws.sensor_pos, ws.refellipsoid, ws.atmosphere_dim, ws.z_tan)
ws.Matrix1ColFromVector(ws.sensor_los, ws.za)
ws.sensor_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
##############################################################################
#
# Absorption
#
##############################################################################
ws.MatrixCreate("Ja")
ws.MatrixCreate("Jp")
ws.StringCreate("info")
# Species (just O3)
#
# Retrieve for a grid rougly matching 16 to 64 km in 2 km steps.
#
ws.VectorCreate("retrieval_grid")
ws.VectorNLogSpace(ws.retrieval_grid, 25, 10000.0, 10.0)
ws.StringSet(ws.info, "O3 rel analytical")
ws.Print(ws.info, 0)
ws.jacobianInit()
ws.jacobianAddAbsSpecies(
    ws.jacobian_quantities,
    ws.jacobian_agenda,
    ws.atmosphere_dim,
    ws.p_grid,
    ws.lat_grid,
    ws.lon_grid,
    ws.retrieval_grid,
    ws.lat_grid,
    ws.lon_grid,
    "O3",
    "rel",
    1,
)
ws.jacobianClose()
# WriteXML( in=jacobian_quantities, filename="Jq_O3_analytical.xml" )
# No scattering
#
ws.cloudboxOff()
#
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.yCalc()
# WriteXML( in=y, filename="y_O3_analytical.xml" )
# WriteXML( in=jacobian, filename="J_O3_analytical.xml" )
ws.Copy(ws.Ja, ws.jacobian)
# Same with perturbations
#
ws.StringSet(ws.info, "O3 rel perturbation")
ws.Print(ws.info, 0)
#
ws.IndexNumberOfAtmosphericPoints(n=ws.ybatch_n, p_grid=ws.retrieval_grid)
ws.NumericCreate("perturbation")
ws.NumericSet(ws.perturbation, 0.01)
ws.jacobianOff()
#
@arts_agenda
def ybatch_calc_agenda(ws):
    ws.vmr_fieldPerturb(
        p_ret_grid=ws.retrieval_grid,
        lat_ret_grid=ws.lat_grid,
        lon_ret_grid=ws.lon_grid,
        species="O3",
        pert_index=ws.ybatch_index,
        pert_size=ws.perturbation,
        pert_mode="relative",
    )
    ws.yCalc()


ws.ybatch_calc_agenda = ybatch_calc_agenda

#
ws.ybatchCalc(ybatch_start=0)
ws.jacobianFromYbatch(pert_size=ws.perturbation)
# WriteXML( in=y, filename="y_O3_perturbation.xml" )
# WriteXML( in=jacobian, filename="J_O3_perturbation.xml" )
ws.Copy(ws.Jp, ws.jacobian)
#
ws.Compare(ws.Ja, ws.Jp, 0.005)
##############################################################################
#
# Temperature, without HSE
#
##############################################################################
# For limb sounding, the analytical expressions do not cover all effects
# related to HSE and refraction and HSE must be "off" here to get consistent
# results.
# Stuff needed for HSE
ws.NumericSet(ws.p_hse, 100000.0)
ws.NumericSet(ws.z_hse_accuracy, 1.0)
ws.VectorSet(ws.lat_true, np.array([15.0]))
ws.VectorSet(ws.lon_true, np.array([123.0]))
#
ws.StringSet(ws.info, "T analytical")
ws.Print(ws.info, 0)
ws.jacobianInit()
ws.jacobianAddTemperature(
    ws.jacobian_quantities,
    ws.jacobian_agenda,
    ws.atmosphere_dim,
    ws.p_grid,
    ws.lat_grid,
    ws.lon_grid,
    ws.retrieval_grid,
    ws.lat_grid,
    ws.lon_grid,
    "off",
)
ws.jacobianClose()
# WriteXML( in=jacobian_quantities, filename="Jq_T_analytical.xml" )
ws.yCalc()
# WriteXML( in=y, filename="y_T_analytical.xml" )
# WriteXML( in=jacobian, filename="J_T_analytical.xml" )
ws.Copy(ws.Ja, ws.jacobian)
# Same with perturbations
#
ws.StringSet(ws.info, "T perturbation")
ws.Print(ws.info, 0)
#
ws.NumericSet(ws.perturbation, 0.1)
ws.jacobianOff()
#
@arts_agenda
def ybatch_calc_agenda(ws):
    ws.AtmFieldPerturb(
        perturbed_field=ws.t_field,
        original_field=ws.t_field,
        p_ret_grid=ws.retrieval_grid,
        lat_ret_grid=ws.lat_grid,
        lon_ret_grid=ws.lon_grid,
        pert_index=ws.ybatch_index,
        pert_size=ws.perturbation,
    )
    ws.yCalc()


ws.ybatch_calc_agenda = ybatch_calc_agenda

#
ws.ybatchCalc(ybatch_start=0)
ws.jacobianFromYbatch(pert_size=ws.perturbation)
#
# WriteXML( in=y, filename="y_T_perturbation.xml" )
# WriteXML( in=jacobian, filename="J_T_perturbation.xml" )
ws.Copy(ws.Jp, ws.jacobian)
# Compare
ws.Compare(ws.Ja, ws.Jp, 0.001)
##############################################################################
#
# Pointing
#
##############################################################################
# Sensor time must be specified here
ws.nrowsGet(ws.nrows, ws.sensor_pos)
ws.VectorNLinSpace(ws.sensor_time, ws.nrows, 0.0, 1.0)
ws.StringSet(ws.info, "Pointing recalc")
ws.Print(ws.info, 0)
ws.jacobianInit()
ws.jacobianAddPointingZa(
    ws.jacobian_quantities,
    ws.jacobian_agenda,
    ws.sensor_pos,
    ws.sensor_time,
    0,
    "recalc",
    0.001,
)
ws.jacobianClose()
ws.yCalc()
ws.Copy(ws.Ja, ws.jacobian)
ws.StringSet(ws.info, "Pointing interp")
ws.Print(ws.info, 0)
ws.jacobianInit()
ws.jacobianAddPointingZa(
    ws.jacobian_quantities,
    ws.jacobian_agenda,
    ws.sensor_pos,
    ws.sensor_time,
    0,
    "interp",
    0.001,
)
ws.jacobianClose()
ws.yCalc()
ws.Copy(ws.Jp, ws.jacobian)
# WriteXML( "ascii", Ja, "J_pointing_recalc.xml" )
# WriteXML( "ascii", Jp, "J_pointing_interp.xml" )
# Compare (Note that the WF is for 1 deg change, corresponding to about
# 60 km change in tangent altitude, and 10 K/deg accuracy is OK)
ws.Compare(ws.Ja, ws.Jp, 10)
##############################################################################
#
# Winds
#
##############################################################################
ws.StringSet(ws.info, "Wind/v, analytical")
ws.Print(ws.info, 0)
#
ws.IndexSet(ws.abs_f_interp_order, 3)
#
ws.jacobianInit()
ws.jacobianAddWind(
    ws.jacobian_quantities,
    ws.jacobian_agenda,
    ws.atmosphere_dim,
    ws.p_grid,
    ws.lat_grid,
    ws.lon_grid,
    ws.retrieval_grid,
    ws.lat_grid,
    ws.lon_grid,
    "v",
)
ws.jacobianClose()
# WriteXML( in=jacobian_quantities, filename="Jq_wind.xml" )
ws.yCalc()
ws.Copy(ws.Ja, ws.jacobian)
# WriteXML( in=y, filename="y_wind_analytical.xml" )
# WriteXML( in=jacobian, filename="J_wind_analytical.xml" )
# Same with perturbations
#
# wind_u_field is [], so we need to reset it to zeros of correct size
ws.Copy(ws.wind_v_field, ws.t_field)
ws.Tensor3Scale(ws.wind_v_field, ws.wind_v_field, 0.0)
#
ws.StringSet(ws.info, "Wind/v, perturbation")
ws.Print(ws.info, 0)
#
ws.NumericSet(ws.perturbation, 0.1)
ws.jacobianOff()
#
@arts_agenda
def ybatch_calc_agenda(ws):
    ws.AtmFieldPerturb(
        perturbed_field=ws.wind_v_field,
        original_field=ws.wind_v_field,
        p_ret_grid=ws.retrieval_grid,
        lat_ret_grid=ws.lat_grid,
        lon_ret_grid=ws.lon_grid,
        pert_index=ws.ybatch_index,
        pert_size=ws.perturbation,
    )
    ws.yCalc()


ws.ybatch_calc_agenda = ybatch_calc_agenda

#
ws.ybatchCalc(ybatch_start=0)
ws.jacobianFromYbatch(pert_size=ws.perturbation)
#
# WriteXML( in=y, filename="y_wind_perturbation.xml" )
# WriteXML( in=jacobian, filename="J_wind_perturbation.xml" )
ws.Copy(ws.Jp, ws.jacobian)
# Compare
ws.Compare(ws.Ja, ws.Jp, 1e-06)
##############################################################################
#
# Just check that remaining weighting functions don't cause any error
#
##############################################################################
ws.StringSet(ws.info, "Others: frequency and baseline")
ws.Print(ws.info, 0)
ws.jacobianInit()
ws.jacobianAddFreqShift(df=50000.0)
ws.jacobianAddFreqStretch(df=50000.0)
ws.jacobianAddPolyfit(
    ws.jacobian_quantities,
    ws.jacobian_agenda,
    ws.sensor_response_pol_grid,
    ws.sensor_response_dlos_grid,
    ws.sensor_pos,
    1,
    0,
    0,
    0,
)
ws.jacobianAddSinefit(
    ws.jacobian_quantities,
    ws.jacobian_agenda,
    ws.sensor_response_pol_grid,
    ws.sensor_response_dlos_grid,
    ws.sensor_pos,
    np.array([2.0e08, 4.0e07]),
    0,
    0,
    0,
)
ws.jacobianClose()
# WriteXML( in=jacobian_quantities, filename="Jq_Other.xml" )
ws.yCalc()
# WriteXML( in=y, filename="y_Other.xml" )
# WriteXML( in=jacobian, filename="J_Other.xml" )
##############################################################################
#
# Transmission weighting functions
#
##############################################################################
ws.StringSet(ws.info, "Transmission analytical")
ws.Print(ws.info, 0)
# (standard) transmission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Transmission)
ws.Copy(ws.iy_transmitter_agenda, ws.iy_transmitter_agenda__UnitUnpolIntensity)
# Calculate with analytical
ws.jacobianInit()
ws.jacobianAddTemperature(
    ws.jacobian_quantities,
    ws.jacobian_agenda,
    ws.atmosphere_dim,
    ws.p_grid,
    ws.lat_grid,
    ws.lon_grid,
    ws.retrieval_grid,
    ws.lat_grid,
    ws.lon_grid,
    "off",
)
ws.jacobianClose()
ws.yCalc()
ws.Copy(ws.Ja, ws.jacobian)
# Same with perturbations
#
ws.StringSet(ws.info, "T perturbation")
ws.Print(ws.info, 0)
#
ws.NumericSet(ws.perturbation, 0.1)
ws.jacobianOff()
#
@arts_agenda
def ybatch_calc_agenda(ws):
    ws.AtmFieldPerturb(
        perturbed_field=ws.t_field,
        original_field=ws.t_field,
        p_ret_grid=ws.retrieval_grid,
        lat_ret_grid=ws.lat_grid,
        lon_ret_grid=ws.lon_grid,
        pert_index=ws.ybatch_index,
        pert_size=ws.perturbation,
    )
    ws.yCalc()


ws.ybatch_calc_agenda = ybatch_calc_agenda

#
ws.ybatchCalc(ybatch_start=0)
ws.jacobianFromYbatch(pert_size=ws.perturbation)
#
ws.Copy(ws.Jp, ws.jacobian)
# Compare
ws.Compare(ws.Ja, ws.Jp, 0.001)
