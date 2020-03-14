# example: limb case using an optimized polar angle grid.

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
# No jacobian calculations
#
ws.jacobianOff()
# Frequency grid
# --------------
# Note: The frequencies must be contained in the gas absorption lookup table.
ws.VectorSet(ws.f_grid, np.array([2.295e11, 2.305e11]))
# Number of Stokes components to be computed
# -------------------------------------------
# IndexSet( stokes_dim, 4 )
# Radiance output units
# ----------------------
ws.StringSet(ws.iy_unit, "RJBT")
# Definition of the atmosphere
# ----------------------------
# Dimensionality of the atmosphere
ws.AtmosphereSet1D()
# Pressure grid
ws.ReadXML(ws.p_grid, "testdata/testdoit_p_grid.xml")
# Definition of species
ws.abs_speciesSet(species=["H2O-PWR98", "O2-PWR93", "N2-SelfContStandardType"])
# Atmospheric profiles
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalc()
# Gas absorption from lookup table
# ---------------------------------
# in case an abs_lookup table need to be calculated, uncomment the block below
# (and replace filename for abs_lookup table to read in ReadXML call).
# atmfields_checkedCalc
# abs_lines_per_speciesSetEmpty
# abs_xsec_agenda_checkedCalc
# abs_lookupSetup
# abs_lookupCalc
# WriteXML( output_file_format="binary", in=abs_lookup, filename="my_gas_abs_lookup.xml" )
ws.ReadXML(ws.abs_lookup, "testdata/testdoit_gas_abs_lookup.xml")
ws.abs_lookupAdapt()
# absorption from LUT
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# Definition of Earth surface
# ----------------------------
ws.MatrixSetConstant(ws.z_surface, 1, 1, 500.0)
# Emission and reflection properties:
# - specular reflecting water surface at ambient atmospheric temperature
# - water reflective indices from Liebe model
# - surface skin temperature interpolated from atmospheric t_field
#
ws.VectorCreate("n_t_grid")


@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
    ws.VectorSetConstant(ws.n_t_grid, 1, ws.surface_skin_t)
    ws.complex_refr_indexWaterLiebe93(
        complex_refr_index=ws.surface_complex_refr_index,
        data_f_grid=ws.f_grid,
        data_T_grid=ws.n_t_grid,
    )
    ws.surfaceFlatRefractiveIndex()


ws.surface_rtprop_agenda = surface_rtprop_agenda

# Definition of sensor position and LOS
# --------------------------------------
ws.VectorCreate("vector_1")
ws.VectorCreate("vector_2")
# Sensor viewing angles
ws.VectorSet(ws.vector_1, np.array([99.7841942]))
# VectorNLinSpace( vector_2, 4, 120, 180 ) #extend with a couple of downviewing los
# Append( vector_1, vector_2 )
ws.nelemGet(ws.nelem, ws.vector_1)
# Sensor altitude from earth surface
ws.VectorSetConstant(ws.vector_2, ws.nelem, 95000.1)
ws.Matrix1ColFromVector(ws.sensor_pos, ws.vector_2)
ws.Matrix1ColFromVector(ws.sensor_los, ws.vector_1)
# SensorOff means that the result of the calculation are the radiances,
# which are not modified by sensor properties
ws.sensorOff()
# Specification of cloud
# -----------------------
# Set the cloudbox limits (pressure units)
# Alternative: cloudboxSetManuallyAltitude (specification in [km])
# ----------------------------------------------------------------
ws.cloudboxSetManually(
    p1=71617.7922264, p2=17111.6808705, lat1=0.0, lat2=0.0, lon1=0.0, lon2=0.0
)
# Initialization and adding scattering elements.
# ----------------------------------------------------------------
ws.ScatSpeciesInit()
# Adding scattering elements.
# ----------------------------------------------------------------
# Here actually both added elements are indentical. however, for testing and for
# demonstration purposed, having 2 elements is better.
ws.ScatElementsPndAndScatAdd(
    scat_data_files=[
        "testdata/scatData/azi-random_f229-231T214-225r100NP-1ar1_5ice.xml",
        "testdata/scatData/azi-random_f229-231T214-225r100NP-1ar1_5ice.xml",
    ],
    pnd_field_files=[
        "testdata/testdoit_pnd_field_1D.xml",
        "testdata/testdoit_pnd_field_1D.xml",
    ],
)
# scat_dataCheck( scat_data = scat_data_raw )
ws.scat_dataCalc()
# scat_dataCheck
ws.pnd_fieldCalcFrompnd_field_raw()
# rescaling back to pnd equivalent to only one scattering element. in a true
# case, this is NOT needed!
ws.Tensor4Scale(ws.pnd_field, ws.pnd_field, 0.5)
# -------------
# DOIT settings
# -------------
# Select interpolation method ('linear' or 'polynomial'):
# ----------------------------------------------------
# For limb calculations is is very important to have a fine resolution
# about 90°.
# If "polynomial" is selected one has to use an optimized grid. Please
# use *doit_za_grid_optCalc* to optimize the grid.
ws.doit_za_interpSet(interp_method="linear")
# As one needs for the RT calculations in limb direction a much finer
# zenith angle grid resolution as for the computation of the scattering
# integral, it is necessary to define two zenith angle grids. For the scattering
# integral equidistant grids are created from N_za_grid, N_aa_grid, which give
# the number of grid points. An optimized grid for the RT calculation needs to
# be given by a file. If no filename is specified (za_grid_opt_file = ""), the
# equidistant grids are used for both, scattering integral and RT calculation.
# This option can only be used for down-looking geometries.
ws.DOAngularGridsSet(
    N_za_grid=19, N_aa_grid=37, za_grid_opt_file="testdata/testdoit_za_grid_opt.xml"
)
# Main agenda for DOIT calculation
# --------------------------------
#
# Input: incoming field on the cloudbox boundary
# Ouput: the scattered field on the cloudbox boundary
@arts_agenda
def doit_mono_agenda(ws):
    # Prepare scattering data for DOIT calculation (Optimized method):
    ws.DoitScatteringDataPrepare()
    ws.Ignore(ws.f_grid)
    # Alternative method (needs less memory):
    # scat_data_monoCalc
    # Perform iterations: 1. scattering integral. 2. RT calculations with
    # fixed scattering integral field, 3. convergence test
    ws.cloudbox_field_monoIterate()
    # Write the radiation field inside the cloudbox:
    # WriteXMLIndexed( in=cloudbox_field_mono, file_index=f_index )


ws.doit_mono_agenda = doit_mono_agenda

# Definitions for methods used in *i_fieldIterate*:
# ----------------------------------------------------
# 1. Scattering integral
# --------------------------
# Calculation of the phase matrix
@arts_agenda
def pha_mat_spt_agenda(ws):
    # Optimized option:
    ws.pha_mat_sptFromDataDOITOpt()
    # Alternative option:
    # pha_mat_sptFromMonoData


ws.pha_mat_spt_agenda = pha_mat_spt_agenda


@arts_agenda
def doit_scat_field_agenda(ws):
    ws.doit_scat_fieldCalcLimb()
    # Alternative: use the same za grids in RT part and scattering integral part
    # doit_scat_fieldCalc


ws.doit_scat_field_agenda = doit_scat_field_agenda

# 2. Radiative transfer with fixed scattering integral term
# ---------------------------------------------------------
@arts_agenda
def doit_rte_agenda(ws):
    # Sequential update for 1D
    ws.cloudbox_fieldUpdateSeq1D(normalize=1, norm_error_threshold=0.05)
    # Alternatives:
    # Plane parallel approximation for determination of propagation path steps
    # cloudbox_fieldUpdateSeq1DPP
    # Without sequential update (not efficient):
    # cloudbox_fieldUpdate1D
    # 3D atmosphere:
    # cloudbox_fieldUpdateSeq3D


ws.doit_rte_agenda = doit_rte_agenda

# Calculate opticle properties of particles and add particle absorption
# and extiction to the gaseous properties to get total extinction and
# absorption:
@arts_agenda
def spt_calc_agenda(ws):
    ws.opt_prop_sptFromMonoData()


ws.spt_calc_agenda = spt_calc_agenda

# 3. Convergence test
# ----------------------
@arts_agenda
def doit_conv_test_agenda(ws):
    # Give limits for all Stokes components in Rayleigh Jeans BT:
    ws.doit_conv_flagAbsBT(epsilon=np.array([0.1, 0.01, 0.01, 0.01]))
    # Alternative: Give limits in radiances
    # doit_conv_flagAbs( doit_conv_flag, doit_iteration_counter, cloudbox_field,
    #                   cloudbox_field_old ){
    #  epsilon = [0.1e-15, 0.1e-18, 0.1e-18, 0.1e-18]
    # }
    # If you want to investigat several iteration fields, for example
    # to investigate the convergence behavior, you can use
    # the following method:
    # DoitWriteIterationFields
    # Print( doit_iteration_counter, 0 )


ws.doit_conv_test_agenda = doit_conv_test_agenda


@arts_agenda
def iy_cloudbox_agenda(ws):
    ws.iyInterpCloudboxField()


ws.iy_cloudbox_agenda = iy_cloudbox_agenda

# End of Main
