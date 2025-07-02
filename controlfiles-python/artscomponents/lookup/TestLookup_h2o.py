import pyarts
import matplotlib.pyplot as plt
import numpy as np


ws = pyarts.Workspace()

ws.water_p_eq_agendaSet()
ws.gas_scattering_agendaSet()
ws.PlanetSet(option="Earth")

ws.abs_speciesSet(species=["H2O-PWR98"])
ws.abs_lines_per_speciesSetEmpty()

ws.abs_lookupSetupWide(p_step=0.5)

ws.VectorNLinSpace(output=ws.f_grid, nelem=5, start=1e9, stop=1e12)

ws.stokes_dim = 1

ws.propmat_clearsky_agendaAuto()
ws.lbl_checkedCalc()
ws.jacobianOff()

ws.abs_lookupCalc()

# ws.abs_lookup.value.savexml("TestLookup_h2oREFERENCE.xml", type="zascii")
ref_lookup = pyarts.arts.GasAbsLookup.fromxml("TestLookup_h2oREFERENCE.xml")
# print(np.max(np.abs(ref_lookup.xsec - ws.abs_lookup.value.xsec)))
assert ref_lookup.xsec.shape == ws.abs_lookup.value.xsec.shape
assert np.allclose(ref_lookup.xsec, ws.abs_lookup.value.xsec, atol=1e-15)

# pyarts.plots.plot_arts_lookup(ref_lookup, opacity=False)
# pyarts.plots.plot_arts_lookup(ws.abs_lookup.value, opacity=False)
# plt.show()
