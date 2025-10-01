import pyarts3 as pyarts
import numpy as np


ws = pyarts.Workspace()

ws.atmospheric_fieldInit(toa=100e3)
ws.surface_fieldEarth()
ws.ray_path_observer_agendaSetGeometric()

ws2 = ws.ray_path_observer_agendaExecute(
  spectral_radiance_observer_position = [300e3, 0, 0],
  spectral_radiance_observer_line_of_sight = [180.0, 0.0]
)

ws2.atmospheric_field["t"] = 123.4

assert ws.atmospheric_field["t"].data == 123.4, "Shared data not working"
