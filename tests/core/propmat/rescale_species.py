import pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

ws.abs_speciesSet(species=["O2-PWR98"])
ws.propmat_clearsky_agenda_checked = True

@pyarts.workspace.arts_agenda(ws=ws, fix=True)
def propmat_clearsky_agenda(ws):
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddPredefined()
print("First agenda")
for m in propmat_clearsky_agenda.methods:
    print(m, end='\n\n')

ws.jacobianOff()
ws.f_grid = np.linspace(30e9, 1e12, 1001)
ws.atm_point.temperature = 200
ws.atm_point.pressure = 1e4
ws.atm_point[ws.abs_species[0]] = 0.21
ws.atm_point.mag = [0., 0, 0.]
ws.rtp_los = [0, 0]

ws.propmat_clearsky_agendaExecute()
orig = ws.propmat_clearsky * 1.0

@pyarts.workspace.arts_agenda(ws=ws)
def propmat_clearsky_agenda(ws):
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddPredefined()
    ws.propmat_clearskyAddScaledSpecies(target="O2-PWR98", scale=0.1)
print("Second agenda")
for m in propmat_clearsky_agenda.methods:
    print(m, end='\n\n')


ws.propmat_clearsky_agendaExecute()
scaled = ws.propmat_clearsky * 1.0

assert np.allclose(scaled[:, 0] / orig[:, 0], 1.1), "Something is wrong"
