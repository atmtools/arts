# DEFINITIONS:  -*-sh-*-
#
# Demonstration of a DOIT scattering calculation initialized from a previous
# solution instead of the clearsky field.
#
# Compared to TestDOIT the pnd_field is slightly varied and
# DOIT is initialized by the solution of TestDOIT.
#
# Author: Jana Mendrok
#

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.Tensor7Create("cloudbox_field_ref")
ws.IndexSet(ws.stokes_dim, 4)
ws.execute_controlfile("artscomponents/doit/doit_setup.arts")
# Perturbing the pnd_field by 1%
ws.Tensor4Scale(ws.pnd_field, ws.pnd_field, 1.01)
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
    # WriteXMLIndexed( in=cloudbox_field_mono, file_index=f_index )
    # Perform iterations: 1. scattering integral. 2. RT calculations with
    # fixed scattering integral field, 3. convergence test
    ws.cloudbox_field_monoIterate()
    # Write the radiation field inside the cloudbox:
    # WriteXMLIndexed( in=cloudbox_field_mono, file_index=f_index )


ws.doit_mono_agenda = doit_mono_agenda

# DOIT calculation
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.scat_data_checkedCalc()
ws.sensor_checkedCalc()
# Initialize Doit variables
ws.DoitInit()
# Calculate incoming radiation field at cloudbox boundary
ws.DoitGetIncoming()
ws.ReadXML(ws.cloudbox_field_ref, "doit_i_fieldREFERENCE_DOIT.xml")
ws.cloudbox_fieldSetFromPrecalc(cloudbox_field_precalc=ws.cloudbox_field_ref)
# Executes doit_mono_agenda for all frequencies
ws.DoitCalc()
# Calculate RT from cloudbox boundary to the sensor
#
ws.StringSet(ws.iy_unit, "RJBT")
ws.yCalc()
ws.WriteXML(ws.output_file_format, ws.y, "", 0)
# ==================stop==========================
# ==================check==========================
ws.VectorCreate("yREFERENCE")
ws.ReadXML(ws.yREFERENCE, "artscomponents/doit/yREFERENCE_DOITprecalcInit.xml")
ws.Compare(ws.y, ws.yREFERENCE, 0.2)
# End of Main
