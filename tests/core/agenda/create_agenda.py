import pyarts3 as pyarts

ws = pyarts.Workspace()


@pyarts.arts_agenda
def create_agenda(ws: pyarts.Workspace):
    ws.abs_bandsFromModelState()


@pyarts.in_parallel(ws=ws)
def run_in_parallel(ws: pyarts.Workspace):
    ws.abs_speciesSet(species=["H2O", "CO2"])
    ws.atm_pointInit()
    ws.measurement_sensorsInit()
    ws.measurement_sensorsAddSimple(freq_grid=[1, 2], pos=[3, 4, 5], los=[6, 7])
    ws.measurement_sensorsAddSimple(freq_grid=[1, 2], pos=[3, 4, 5], los=[6, 7])


assert ws.abs_species == ["H2O", "CO2"]
assert len(ws.measurement_sensors) == 2

assert len(run_in_parallel.par_tasks(ws)) == 3
assert len(create_agenda.par_tasks(ws)) == 1
