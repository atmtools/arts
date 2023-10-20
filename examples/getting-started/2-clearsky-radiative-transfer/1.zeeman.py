import pyarts
import os
import numpy as np

ws = pyarts.workspace.Workspace()

# %% Calculations

ws.jacobianOff()
ws.iy_unit = "PlanckBT"
ws.ppath_lstep = 5e2
ws.ppath_lmax = 1e3
ws.ppath_lraytrace = 1e3
ws.rt_integration_option = "default"
ws.rte_alonglos_v = 0.0
ws.nlteOff()

# %% Sensor

ws.f_grid = np.linspace(-4e6, 4e6, 1001) + 118750348044.712
ws.rte_pos = [300e3, 0, 0]
ws.rte_los = [180, 0]

# %% Species and line absorption

ws.abs_speciesSet(species=[f"O2-Z-66-118e9-120e9"])
ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename = "lines/")
ws.Wigner6Init()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propmat_clearsky_agendaAuto()

# %% Grids and planet

ws.PlanetSet(option = "Earth")
ws.atm_fieldInit(toa=100e3)
ws.atm_fieldAddGriddedData(key=pyarts.arts.String("t"),
                           data=pyarts.arts.GriddedField3.fromxml("planets/Earth/afgl/tropical/t.xml"))
ws.atm_fieldAddGriddedData(key=pyarts.arts.String("p"),
                           data=pyarts.arts.GriddedField3.fromxml("planets/Earth/afgl/tropical/p.xml"))
ws.atm_fieldAddGriddedData(key=ws.abs_species[0],
                           data=pyarts.arts.GriddedField3.fromxml("planets/Earth/afgl/tropical/O2.xml"))
ws.atm_fieldIGRF(time="2000-03-11 14:39:37")

# %% Checks

ws.lbl_checked = 1

# %% Core calculations
ws.ppathGeometric()
ws.ppvar_atmFromPath()
ws.ppvar_fFromPath()
ws.ppvar_propmatCalc()
ws.ppvar_srcFromPropmat()
ws.MatrixPlanck(ws.iy, f=ws.f_grid, t=293)
ws.background_radFromMatrix(iy_mat=ws.iy)
ws.ppvar_tramatCalc()
ws.ppvar_cumtramatForward()
ws.ppvar_radCalcEmission()


ws.iyCopyPath()
ws.Touch(ws.diy_dx)
ws.iyUnitConversion(iy_agenda_call1=1)

