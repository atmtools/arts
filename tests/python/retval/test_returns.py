import pyarts3 as pyarts
import numpy as np


ws = pyarts.Workspace()

ws.atm_fieldInit(toa=100e3)
ws.surf_fieldEarth()
ws.ray_path_observer_agendaSetGeometric()

ws2 = ws.ray_path_observer_agendaExecute(
  obs_pos = [300e3, 0, 0],
  obs_los = [180.0, 0.0]
)

ws2.atm_field["t"] = 123.4

assert ws.atm_field["t"].data == 123.4, "Shared data not working"
