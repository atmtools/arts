# DEFINITIONS:  -*-sh-*-
#
# To be written ...
#

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
# on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# Blackbody surface
ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Blackbody_SurfTFromt_field)
ws.VectorSet(ws.surface_scalar_reflectivity, np.array([0.0]))
# Standard ppath calculations
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# Radiative transfer agendas
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
ws.Copy(ws.iy_cloudbox_agenda, ws.iy_cloudbox_agenda__QuadInterpField1D)
# Absorption species
ws.abs_speciesSet(species=["N2-SelfContStandardType", "O2-PWR93", "H2O-PWR98"])
# No line data needed here
ws.abs_lines_per_speciesSetEmpty()
# Dimensionality of the atmosphere
ws.AtmosphereSet1D()
# Brigtness temperatures used
ws.StringSet(ws.iy_unit, "PlanckBT")
# Various things not used
ws.ArrayOfStringSet(ws.iy_aux_vars, [])
ws.jacobianOff()
# Read data created by setup_test.m
ws.ReadXML(ws.p_grid, "testdata/p_grid.xml")
ws.ReadXML(ws.t_field, "testdata/t_field.xml")
ws.ReadXML(ws.z_field, "testdata/z_field.xml")
ws.ReadXML(ws.vmr_field, "testdata/vmr_field.xml")
ws.ReadXML(ws.particle_bulkprop_field, "testdata/particle_bulkprop_field")
ws.ReadXML(ws.particle_bulkprop_names, "testdata/particle_bulkprop_names")
ws.ReadXML(ws.scat_data_raw, "testdata/scat_data.xml")
ws.ReadXML(ws.scat_meta, "testdata/scat_meta.xml")
# Define hydrometeors
#
ws.StringCreate("species_id_string")
#
# Scat species 0
ws.StringSet(ws.species_id_string, "RWC")
ws.ArrayOfStringSet(ws.pnd_agenda_input_names, ["RWC"])


@arts_agenda
def pnd_agenda_array(ws):
    ws.ScatSpeciesSizeMassInfo(species_index=ws.agenda_array_index, x_unit="dveq")
    ws.Copy(ws.psd_size_grid, ws.scat_species_x)
    ws.Copy(ws.pnd_size_grid, ws.scat_species_x)
    ws.psdWangEtAl16(t_min=273.0, t_max=999.0)
    ws.pndFromPsdBasic()


ws.Append(ws.pnd_agenda_array, pnd_agenda_array)

ws.Append(ws.scat_species, ws.species_id_string)
ws.Append(ws.pnd_agenda_array_input_names, ws.pnd_agenda_input_names)
#
# Scat species 1
ws.StringSet(ws.species_id_string, "IWC")
ws.ArrayOfStringSet(ws.pnd_agenda_input_names, ["IWC"])


@arts_agenda
def pnd_agenda_array(ws):
    ws.ScatSpeciesSizeMassInfo(
        species_index=ws.agenda_array_index, x_unit="dveq", x_fit_start=0.0001
    )
    ws.Copy(ws.psd_size_grid, ws.scat_species_x)
    ws.Copy(ws.pnd_size_grid, ws.scat_species_x)
    ws.psdMcFarquaharHeymsfield97(t_min=10.0, t_max=273.0, t_min_psd=210.0)
    ws.pndFromPsdBasic()


ws.Append(ws.pnd_agenda_array, pnd_agenda_array)

ws.Append(ws.scat_species, ws.species_id_string)
ws.Append(ws.pnd_agenda_array_input_names, ws.pnd_agenda_input_names)
# Special settings for the scattering solvers
#
# Angular grids for DOIT and DISORT
ws.DOAngularGridsSet(N_za_grid=38, N_aa_grid=37)
#
# DOIT stuff
@arts_agenda
def pha_mat_spt_agenda(ws):
    ws.pha_mat_sptFromDataDOITOpt()


ws.pha_mat_spt_agenda = pha_mat_spt_agenda


@arts_agenda
def spt_calc_agenda(ws):
    ws.opt_prop_sptFromMonoData()


ws.spt_calc_agenda = spt_calc_agenda


@arts_agenda
def doit_scat_field_agenda(ws):
    ws.doit_scat_fieldCalc()


ws.doit_scat_field_agenda = doit_scat_field_agenda


@arts_agenda
def doit_mono_agenda(ws):
    ws.DoitScatteringDataPrepare()
    ws.Ignore(ws.f_grid)
    ws.cloudbox_field_monoIterate()


ws.doit_mono_agenda = doit_mono_agenda


@arts_agenda
def doit_rte_agenda(ws):
    ws.cloudbox_fieldUpdateSeq1D(normalize=1, norm_error_threshold=0.05)


ws.doit_rte_agenda = doit_rte_agenda

ws.doit_za_interpSet(interp_method="linear")


@arts_agenda
def doit_conv_test_agenda(ws):
    ws.doit_conv_flagAbsBT(epsilon=np.array([0.1]))


ws.doit_conv_test_agenda = doit_conv_test_agenda

#
# RT4 creates own grirs, so we need copies of the ones created above
ws.VectorCreate("za_grid_copy")
ws.VectorCreate("aa_grid_copy")
#
# Hybrid requires that ppath_lmax is not too high
ws.NumericSet(ws.ppath_lmax, 100.0)
ws.AgendaCreate("iy_hybrid_agenda")


@arts_agenda
def iy_hybrid_agenda(ws):
    ws.Ignore(ws.iy_id)
    ws.ppathCalc(cloudbox_on=0)
    ws.iyHybrid()
    # The line below is just temporary
    ws.Touch(ws.iy_aux)


ws.iy_hybrid_agenda = iy_hybrid_agenda

# Versions of y for various calculations
ws.VectorCreate("y_doit")
ws.VectorCreate("y_disort")
ws.VectorCreate("y_hybrid")
ws.VectorCreate("y_rt4")
#
ws.StringCreate("message")
# Perform some basic checks
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
# Intitial settings for tests
ws.IndexSet(ws.stokes_dim, 1)
# Scattering data tailored to these frequencies, so don't change!
ws.VectorSet(ws.f_grid, np.array([3.15e10, 1.65e11, 6.66e11]))
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.MatrixSet(ws.sensor_pos, np.array([[20000.0], [20000.0], [10000.0], [5000.0]]))
ws.MatrixSet(ws.sensor_los, np.array([[180.0], [130.0], [160.0], [20.0]]))
# Some stuff that depends on the settings above
ws.sensorOff()
ws.atmgeom_checkedCalc()
ws.sensor_checkedCalc()
ws.scat_dataCalc()
# We need here a higher sca_mat_threshold. Otherwise there is an error for RWC
# and 668 GHz
ws.scat_data_checkedCalc(sca_mat_threshold=0.25)
#
ws.VectorExtractFromMatrix(ws.rtp_pos, ws.z_surface, 0, "row")
ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
# Test 1: Cloubox with no scattering
# ---------------------------------------------------------------------
ws.cloudboxSetFullAtm()
ws.pnd_fieldZero()
ws.cloudbox_checkedCalc()
#
ws.StringSet(ws.message, "No scattering")
ws.Print(ws.message, 0)
#
ws.execute_controlfile("run_doit.arts")
ws.Print(ws.y_doit, 0)
ws.execute_controlfile("run_rt4.arts")
ws.Print(ws.y_rt4, 0)
ws.execute_controlfile("run_disort.arts")
ws.Print(ws.y_disort, 0)
ws.execute_controlfile("run_hybrid.arts")
ws.Print(ws.y_hybrid, 0)
#
ws.Compare(ws.y_doit, ws.y_hybrid, 0.1, "Zero particles, DOIT")
ws.Compare(ws.y_rt4, ws.y_hybrid, 0.1, "Zero particles, RT4")
ws.Compare(ws.y_disort, ws.y_hybrid, 0.1, "Zero particles, DISORT")
# Test 2: With nominal RWC/IWC and surface at p_grid[0]
# ---------------------------------------------------------------------
ws.pnd_fieldCalcFromParticleBulkProps()
ws.cloudbox_checkedCalc()
#
ws.StringSet(ws.message, "Nominal case")
ws.Print(ws.message, 0)
#
ws.execute_controlfile("run_doit.arts")
ws.Print(ws.y_doit, 0)
ws.execute_controlfile("run_rt4.arts")
ws.Print(ws.y_rt4, 0)
ws.execute_controlfile("run_disort.arts")
ws.Print(ws.y_disort, 0)
ws.execute_controlfile("run_hybrid.arts")
ws.Print(ws.y_hybrid, 0)
#
ws.Compare(ws.y_doit, ws.y_disort, 0.4, "Test2, DOIT vs. DISORT")
ws.Compare(ws.y_rt4, ws.y_disort, 0.8, "Test2, RT4 vs. DISORT")
ws.Compare(ws.y_hybrid, ws.y_disort, 0.4, "Test2, Hybrid vs. DISORT")
# Test 3: As 2 but with RWC/IWC increased with a factor of 3
# ---------------------------------------------------------------------
ws.Tensor4Scale(ws.particle_bulkprop_field, ws.particle_bulkprop_field, 3.0)
ws.pnd_fieldCalcFromParticleBulkProps()
#
ws.StringSet(ws.message, "Increased RWC/IWC")
ws.Print(ws.message, 0)
#
ws.execute_controlfile("run_doit.arts")
ws.Print(ws.y_doit, 0)
ws.execute_controlfile("run_rt4.arts")
ws.Print(ws.y_rt4, 0)
ws.execute_controlfile("run_disort.arts")
ws.Print(ws.y_disort, 0)
ws.execute_controlfile("run_hybrid.arts")
ws.Print(ws.y_hybrid, 0)
#
ws.Compare(ws.y_doit, ws.y_disort, 1.0, "Test3, DOIT vs. DISORT")
ws.Compare(ws.y_rt4, ws.y_disort, 0.8, "Test3, RT4 vs. DISORT")
ws.Compare(ws.y_hybrid, ws.y_disort, 0.4, "Test3, Hybrid vs. DISORT")
# Test 4: As 2 but with surface moved upwards
# ---------------------------------------------------------------------
ws.Tensor4Scale(ws.particle_bulkprop_field, ws.particle_bulkprop_field, 0.3333333)
ws.pnd_fieldCalcFromParticleBulkProps()
ws.MatrixAddScalar(ws.z_surface, ws.z_surface, 3000.0)
ws.VectorExtractFromMatrix(ws.rtp_pos, ws.z_surface, 0, "row")
ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
##
ws.StringSet(ws.message, "Surface moved upwards")
ws.Print(ws.message, 0)
#
ws.execute_controlfile("run_rt4.arts")
ws.Print(ws.y_rt4, 0)
ws.execute_controlfile("run_disort.arts")
ws.Print(ws.y_disort, 0)
ws.execute_controlfile("run_hybrid.arts")
ws.Print(ws.y_hybrid, 0)
#
ws.Compare(ws.y_rt4, ws.y_disort, 0.8, "Test4, RT4 vs. DISORT")
ws.Compare(ws.y_hybrid, ws.y_disort, 0.4, "Test4, Hybrid vs. DISORT")
