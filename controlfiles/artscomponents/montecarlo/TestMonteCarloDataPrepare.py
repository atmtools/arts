# DEFINITIONS:  -*-sh-*-
# This control file prepares a lot of atmospheric field data for the ARTS-MC examples
# simpleMCGeneral.arts, and simpleMCGeneralGaussian.arts, simpleMC.arts

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
ws.output_file_formatSetBinary()
ws.VectorSet(ws.f_grid, np.array([2.3e11]))
ws.WriteXML(ws.output_file_format, ws.f_grid)
ws.IndexSet(ws.f_index, 0)
ws.IndexSet(ws.stokes_dim, 4)
ws.AtmosphereSet3D()
ws.ReadXML(ws.p_grid, "p_grid.xml")
ws.ReadXML(ws.lat_grid, "lat_grid.xml")
ws.ReadXML(ws.lon_grid, "lon_grid.xml")
ws.abs_speciesSet(species=["O2-PWR93", "N2-SelfContStandardType", "H2O-PWR98"])
ws.abs_lines_per_speciesSetEmpty()
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalcExpand1D()
ws.WriteXML(ws.output_file_format, ws.t_field)
ws.WriteXML(ws.output_file_format, ws.z_field)
ws.WriteXML(ws.output_file_format, ws.vmr_field)
ws.nelemGet(ws.nrows, ws.lat_grid)
ws.nelemGet(ws.ncols, ws.lon_grid)
ws.MatrixSetConstant(ws.z_surface, ws.nrows, ws.ncols, 500.0)
ws.WriteXML(ws.output_file_format, ws.z_surface)
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.atmfields_checkedCalc()
ws.jacobianOff()
ws.abs_lookupSetup()
ws.abs_lookupCalc()
ws.WriteXML(ws.output_file_format, ws.abs_lookup)
ws.cloudboxSetManually(
    p1=21617.7922264, p2=17111.6808705, lat1=-1.9, lat2=1.9, lon1=-1.9, lon2=1.9
)
ws.Print(ws.cloudbox_limits, 1)
ws.ScatSpeciesInit()
ws.ScatElementsPndAndScatAdd(
    scat_data_files=[
        "testdata/scatData/azi-random_f229-231T214-225r100NP-1ar1_5ice.xml"
    ],
    pnd_field_files=[""],
)
ws.ReadXML(ws.pnd_field_raw, "pnd_field_raw.xml")
ws.pnd_fieldCalcFrompnd_field_raw()
ws.WriteXML(ws.output_file_format, ws.pnd_field)
ws.scat_dataCalc()
# For some reason SingleScatteringData binary files can't be loaded
ws.output_file_formatSetAscii()
ws.WriteXML(ws.output_file_format, ws.scat_data)
ws.WriteXML(ws.output_file_format, ws.cloudbox_limits)
