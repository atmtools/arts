# DEFINITIONS:  -*-sh-*-
# ARTS setup file for DOIT calculations.
#
# This is a control file specifies various Agendas necessary for DOIT
# calculations, and also provides some alternatives.

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Main agenda for DOIT calculation
# ----------------------------------------------------------------------
#
# Input: incoming field on the cloudbox boundary
# Ouput: the scattered field on the cloudbox boundary
@arts_agenda
def doit_mono_agenda(ws):
    # Prepare scattering data for DOIT calculation (Optimized method):
    ws.DoitScatteringDataPrepare()
    ws.Ignore(ws.f_grid)
    # Alternative method:
    # no optimization of scattering angle grids (needs less memory):
    # scat_data_monoCalc
    # Perform iterations: 1. scattering integral. 2. RT calculations with
    # fixed scattering integral field, 3. convergence test
    ws.cloudbox_field_monoIterate()
    # Write the radiation field inside the cloudbox:
    # WriteXML( output_file_format, cloudbox_field )


ws.doit_mono_agenda = doit_mono_agenda

# Definitions for methods used in *i_fieldIterate*:
# --------------------------------------------------
# 1. Scattering integral
# ----------------------
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
    ws.cloudbox_fieldUpdateSeq1D()
    # Alternatives:
    # Plane parallel approximation for determination of propagation path steps:
    # cloudbox_fieldUpdateSeq1DPP
    # Without sequential update (not efficient):
    # cloudbox_fieldUpdate1D
    # 3D atmosphere:
    # cloudbox_fieldUpdateSeq3D
    # Sequential update for 1D with normalizing the scattered field:
    # cloudbox_fieldUpdateSeq1D( normalize=1 )


ws.doit_rte_agenda = doit_rte_agenda

# Calculate opticle properties of particles and add particle absorption
# and extinction to the gaseous properties to get total extinction and
# absorption:
@arts_agenda
def spt_calc_agenda(ws):
    ws.opt_prop_sptFromMonoData()


ws.spt_calc_agenda = spt_calc_agenda

# 3. Convergence test
# --------------------
@arts_agenda
def doit_conv_test_agenda(ws):
    # Give limits for all Stokes components in BT:
    ws.doit_conv_flagAbsBT(epsilon=np.array([0.1]))
    # for stokes dim 1
    # Alternative: Give limits in radiances
    # doit_conv_flagAbs( epsilon=[0.1e-15, 0.1e-18, 0.1e-18, 0.1e-18] )
    #
    # If you want to investigate several iteration fields, for example
    # to investigate the convergence behavior, you can use
    # the following method:
    # DoitWriteIterationFields( iterations=[2, 4])
    # }


ws.doit_conv_test_agenda = doit_conv_test_agenda

ws.Copy(ws.iy_cloudbox_agenda, ws.iy_cloudbox_agenda__LinInterpField)
# End main-----------------------------------------------------------------
