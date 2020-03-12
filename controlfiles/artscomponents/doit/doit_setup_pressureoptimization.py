import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Initialization and adding scattering elements.
# ----------------------------------------------------------------
ws.ScatSpeciesInit()
# Adding scattering elements.
# ----------------------------------------------------------------
ws.ScatElementsPndAndScatAdd(
    scat_data_files=["scattering/H2O_ice/MieSphere_R8.00000e+02um.xml"],
    pnd_field_files=["testdata/testdoit_pnd_field_1D.xml"],
)
ws.scat_dataCheck(scat_data=ws.scat_data_raw)
ws.scat_dataCalc()
ws.scat_dataCheck()
ws.pnd_fieldCalcFrompnd_field_raw()
# here we want a really thin cloud
ws.Tensor4Scale(ws.pnd_field, ws.pnd_field, 0.01)
# Consistency checks needed for running merge
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.ScatSpeciesMerge()
ws.cloudbox_checkedCalc()
ws.scat_data_checkedCalc()
# Main agenda for DOIT calculation
# --------------------------------
@arts_agenda
def doit_mono_agenda(ws):
    # Prepare scattering data for DOIT calculation (Optimized method):
    ws.DoitScatteringDataPrepare()
    ws.Ignore(ws.f_grid)
    # Switch on pressure grid optimization
    ws.OptimizeDoitPressureGrid()
    # Alternative method (needs less memory):
    # scat_data_monoCalc
    # Perform iterations: 1. scattering integral. 2. RT calculations with
    # fixed scattering integral field, 3. convergence test, use acceleration
    ws.cloudbox_field_monoIterate(accelerated=1)
    # Interpolate cloudbox_field_mono back to original size
    ws.cloudbox_field_monoOptimizeReverse()
    # Write the radiation field inside the cloudbox:
    # WriteXMLIndexed( in=cloudbox_field_mono, file_index=f_index )


ws.doit_mono_agenda = doit_mono_agenda

# Convergence test
# ----------------------
@arts_agenda
def doit_conv_test_agenda(ws):
    # Give limits for all Stokes components in Rayleigh Jeans BT:
    ws.doit_conv_flagAbsBT(epsilon=np.array([0.1]))
    # Alternative: Give limits in radiances
    # doit_conv_flagAbs( doit_conv_flag, doit_iteration_counter, cloudbox_field,
    #                   cloudbox_field_old ){
    #  epsilon = [0.1e-15, 0.1e-18, 0.1e-18, 0.1e-18]
    # }
    # If you want to investigat several iteration fields, for example
    # to investigate the convergence behavior, you can use
    # the following method:
    # DoitWriteIterationFields
    ws.Print(ws.doit_iteration_counter, 0)


ws.doit_conv_test_agenda = doit_conv_test_agenda


@arts_agenda
def iy_cloudbox_agenda(ws):
    ws.iyInterpCloudboxField()


ws.iy_cloudbox_agenda = iy_cloudbox_agenda

# End of Main
