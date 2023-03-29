import pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

ws.abs_speciesSet(species=["O2-PWR98"])
ws.abs_lines_per_speciesSetEmpty()
ws.propmat_clearsky_agendaAuto()

ws.jacobianOff()
ws.f_grid = np.linspace(30e9, 1e12, 1001)
ws.stokes_dim = 1
ws.rtp_temperature = 200
ws.rtp_pressure = 1e4
ws.rtp_vmr = [0.21]
ws.rtp_mag = [0, 0, 0]
ws.rtp_los = [0, 0]
ws.Touch(ws.rtp_nlte)

ws.propmat_clearsky_agenda.value.execute(ws)
orig = ws.propmat_clearsky.value.data.value.flatten() * 1.0

@pyarts.workspace.arts_agenda(ws=ws, set_agenda=True)
def propmat_clearsky_agenda(ws):
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddPredefined()
    ws.propmat_clearskyAddScaledSpecies(target="O2-PWR98", scale=0.1)


ws.propmat_clearsky_agenda.value.execute(ws)
scaled = ws.propmat_clearsky.value.data.value.flatten() * 1.0

assert np.allclose(scaled / orig, 1.1), "Something is wrong"
