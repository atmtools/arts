import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
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
# Set first guess field
ws.cloudbox_fieldSetClearsky()
# Executes doit_mono_agenda for all frequencies
ws.DoitCalc()
# Writing out the radition field solution from inside cloudbox for
# initialization of perturbed pnd_field calculation from non-clearsky first
# guess
# WriteXML( in=cloudbox_field, filename="cloudbox_fieldREFERENCE_DOIT.xml" )
# WriteXML( in=cloudbox_field )
# Calculate RT from cloudbox boundary to the sensor
ws.yCalc()
# End of Main
