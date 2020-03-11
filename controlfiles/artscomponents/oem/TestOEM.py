# DEFINITIONS:  -*-sh-*-
#
# Author: Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Agendas to use
#
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# Basic settings
#
ws.AtmosphereSet1D()
ws.IndexSet(ws.stokes_dim, 1)
# Variables to define frequency grid and
#
ws.NumericCreate("f0")
ws.NumericCreate("df_start")
ws.NumericCreate("df_end")
ws.VectorCreate("fgrid1")
ws.VectorCreate("fgrid2")
ws.IndexCreate("nf")
#
ws.NumericSet(ws.f0, 110836000000.0)
#
# Create broad, coarse frequency grid
ws.NumericSet(ws.df_start, -300000000.0)
ws.NumericSet(ws.df_end, 300000000.0)
ws.IndexSet(ws.nf, 601)
ws.VectorNLinSpace(ws.fgrid1, ws.nf, ws.df_start, ws.df_end)
#
# Create a narrow, fine frequency grid
ws.NumericSet(ws.df_start, -10000000.0)
ws.NumericSet(ws.df_end, 10000000.0)
ws.IndexSet(ws.nf, 401)
ws.VectorNLinSpace(ws.fgrid2, ws.nf, ws.df_start, ws.df_end)
#
# Create final f_grid
ws.VectorInsertGridPoints(ws.f_grid, ws.fgrid1, ws.fgrid2)
ws.VectorAddScalar(ws.f_grid, ws.f_grid, ws.f0)
# Define spectrometer
#
ws.NumericCreate("f_resolution")
ws.NumericCreate("f_start")
ws.NumericCreate("f_end")
#
ws.NumericSet(ws.f_start, -280000000.0)
ws.NumericSet(ws.f_end, 280000000.0)
ws.NumericSet(ws.f_resolution, 200000.0)
#
ws.NumericAdd(ws.f_start, ws.f_start, ws.f0)
ws.NumericAdd(ws.f_end, ws.f_end, ws.f0)
# Pressure grid
#
ws.IndexCreate("np")
ws.VectorCreate("p_ret_grid")
#
ws.IndexSet(ws.np, 81)
#
ws.VectorNLogSpace(ws.p_grid, 361, 50000.0, 0.1)
ws.VectorNLogSpace(ws.p_ret_grid, ws.np, 50000.0, 0.1)
# Spectroscopy
#
ws.abs_speciesSet(species=["O3"])
#
ws.ReadARTSCAT(ws.abs_lines, "testdata/ozone_line.xml")
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_lines_per_speciesSetNormalization(option="VVH")
ws.abs_lines_per_speciesSetCutoff(option="ByLine", value=750000000000.0)
# Atmosphere (a priori)
#
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalc()
#
ws.MatrixSetConstant(ws.z_surface, 1, 1, 10000.0)
#
ws.VectorSet(ws.lat_true, np.array([10.0]))
ws.VectorSet(ws.lon_true, np.array([123.0]))
# Apply HSE
#
ws.NumericSet(ws.p_hse, 10000.0)
ws.NumericSet(ws.z_hse_accuracy, 0.5)
#
ws.atmfields_checkedCalc()
#
ws.z_fieldFromHSE()
# Sensor pos/los/time
#
ws.MatrixSetConstant(ws.sensor_pos, 1, 1, 15000.0)
ws.MatrixSetConstant(ws.sensor_los, 1, 1, 60.0)
#
ws.VectorSetConstant(ws.sensor_time, 1, 0.0)
# True f_backend
#
ws.VectorLinSpace(ws.f_backend, ws.f_start, ws.f_end, ws.f_resolution)
# Assume a Gaussian channel response
#
ws.VectorCreate("fwhm")
ws.VectorSetConstant(ws.fwhm, 1, ws.f_resolution)
ws.backend_channel_responseGaussian(fwhm=ws.fwhm)
# With a frequency shift retrieval, we must use sensor_response_agenda
#
ws.FlagOn(ws.sensor_norm)
#
@arts_agenda
def sensor_response_agenda(ws):
    ws.AntennaOff()
    ws.sensor_responseInit()
    # Among responses we only include a backend
    ws.sensor_responseBackend()


ws.sensor_response_agenda = sensor_response_agenda

#
ws.AgendaExecute(ws.sensor_response_agenda)
# RT
#
ws.NumericSet(ws.ppath_lmax, -1.0)
ws.StringSet(ws.iy_unit, "RJBT")
# Deactive parts not used (jacobian activated later)
#
ws.jacobianOff()
ws.cloudboxOff()
# Perform tests
#
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.lbl_checkedCalc()
# Simulate "measurement vector"
#
ws.yCalc()
#
# Start on retrieval specific part
#
# Some vaiables
#
ws.VectorCreate("vars")
ws.SparseCreate("sparse_block")
ws.MatrixCreate("dense_block")
# Start definition of retrieval quantities
#
ws.retrievalDefInit()
#
ws.nelemGet(ws.nelem, ws.p_ret_grid)
ws.nelemGet(ws.nf, ws.sensor_response_f)
# Add ozone as retrieval quantity
#
ws.retrievalAddAbsSpecies(
    species="O3", unit="vmr", g1=ws.p_ret_grid, g2=ws.lat_grid, g3=ws.lon_grid
)
#
ws.VectorSetConstant(ws.vars, ws.nelem, 1e-12)
ws.DiagonalMatrix(ws.sparse_block, ws.vars)
ws.covmat_sxAddBlock(block=ws.sparse_block)
# Add a frquency shift retrieval
#
ws.retrievalAddFreqShift(df=50000.0)
#
ws.VectorSetConstant(ws.vars, 1, 10000000000.0)
ws.DiagonalMatrix(ws.sparse_block, ws.vars)
ws.covmat_sxAddBlock(block=ws.sparse_block)
# Add a baseline fit
#
ws.retrievalAddPolyfit(poly_order=0)
#
ws.VectorSetConstant(ws.vars, 1, 0.5)
ws.DiagonalMatrix(ws.sparse_block, ws.vars)
ws.covmat_sxAddBlock(block=ws.sparse_block)
# Define Se and its invers
#
ws.VectorSetConstant(ws.vars, ws.nf, 0.01)
ws.DiagonalMatrix(ws.sparse_block, ws.vars)
ws.covmat_seAddBlock(block=ws.sparse_block)
#
ws.VectorSetConstant(ws.vars, ws.nf, 100.0)
ws.DiagonalMatrix(ws.dense_block, ws.vars)
ws.covmat_seAddInverseBlock(block=ws.dense_block)
# End definition of retrieval quantities
#
ws.retrievalDefClose()
# x, jacobian and yf must be initialised (or pre-calculated as shown below)
#
ws.VectorSet(ws.x, [])
ws.VectorSet(ws.yf, [])
ws.MatrixSet(ws.jacobian, [])
# Or to pre-set x, jacobian and yf
#
# Copy( x, xa )
# MatrixSet( jacobian, [] )
# AgendaExecute( inversion_iterate_agenda )
# Iteration agenda
#
@arts_agenda
def inversion_iterate_agenda(ws):
    ws.Ignore(ws.inversion_iteration_counter)
    # Map x to ARTS' variables
    ws.x2artsAtmAndSurf()
    ws.x2artsSensor()
    # No need to call this WSM if no sensor variables retrieved
    # To be safe, rerun some checks
    ws.atmfields_checkedCalc()
    ws.atmgeom_checkedCalc()
    # Calculate yf and Jacobian matching x.
    ws.yCalc(y=ws.yf)
    # Add baseline term (no need to call this WSM if no sensor variables retrieved)
    ws.VectorAddVector(ws.yf, ws.yf, ws.y_baseline)
    # This method takes cares of some "fixes" that are needed to get the Jacobian
    # right for iterative solutions. No need to call this WSM for linear inversions.
    ws.jacobianAdjustAndTransform()


ws.inversion_iterate_agenda = inversion_iterate_agenda

# Let a priori be off with 0.5 ppm
#
ws.Tensor4AddScalar(ws.vmr_field, ws.vmr_field, 5e-07)
# Add a baseline
#
ws.VectorAddScalar(ws.y, ws.y, 1.0)
# Introduce a frequency error
#
ws.VectorAddScalar(ws.f_backend, ws.f_backend, -150000.0)
# Calculate sensor_reponse (this time with assumed f_backend)
#
ws.AgendaExecute(ws.sensor_response_agenda)
# Create xa
#
ws.xaStandard()
# Run OEM
ws.OEM(
    method="gn",
    max_iter=5,
    display_progress=1,
    stop_dx=0.1,
    lm_ga_settings=np.array([10.0, 2.0, 2.0, 100.0, 1.0, 99.0]),
)
#
ws.Print(ws.oem_errors, 0)
ws.Print(ws.x, 0)
# Compute averaging kernel matrix
#
ws.avkCalc()
# Compute smoothing error covariance matrix
#
ws.covmat_ssCalc()
# Compute observation system error covariance matrix
#
ws.covmat_soCalc()
# Extract observation errors
#
ws.retrievalErrorsExtract()
# WriteXML( "ascii", f_backend, "f.xml" )
# WriteXML( "ascii", y, "y.xml" )
