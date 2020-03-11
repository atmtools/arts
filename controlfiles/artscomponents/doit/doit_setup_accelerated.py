# setup additions/modifications for DOIT accelerated mode

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Main agenda for DOIT calculation
# --------------------------------
@arts_agenda
def doit_mono_agenda(ws):
    # Prepare scattering data for DOIT calculation (Optimized method):
    ws.DoitScatteringDataPrepare()
    ws.Ignore(ws.f_grid)
    # Alternative method (needs less memory):
    # scat_data_monoCalc
    # Perform iterations: 1. scattering integral. 2. RT calculations with
    # fixed scattering integral field, 3. convergence test
    # Accelerate using NG-Acceleration. To accelerate only the first Stokes component,
    # write accelerated = 1. To accelerate all 4 components, write accelerated = 4.
    # default is set to 0 (no acceleration)
    ws.cloudbox_field_monoIterate(accelerated=4)
    # Write the radiation field inside the cloudbox:
    # WriteXMLIndexed( in=cloudbox_field_mono, file_index=f_index )


ws.doit_mono_agenda = doit_mono_agenda

# Convergence test
# ----------------------
@arts_agenda
def doit_conv_test_agenda(ws):
    # Give limits for all Stokes components in Rayleigh Jeans BT:
    ws.doit_conv_flagAbsBT(epsilon=array([0.001, 0.01, 0.01, 0.01]))
    # Alternative: Give limits in radiances
    # doit_conv_flagAbs( doit_conv_flag, doit_iteration_counter, cloudbox_field,
    #                   cloudbox_field_old ){
    #  epsilon = [0.1e-15, 0.1e-18, 0.1e-18, 0.1e-18]
    # }
    # If you want to investigate several iteration fields, for example
    # to investigate the convergence behavior, you can use
    # the following method:
    # DoitWriteIterationFields
    ws.Print(ws.doit_iteration_counter, 0)


ws.doit_conv_test_agenda = doit_conv_test_agenda

# we want a really thick cloud here
ws.Tensor4Scale(ws.pnd_field, ws.pnd_field, 80.0)
# End of Main
