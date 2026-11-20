"""Operators made of simple structs survive saving and loading as XML."""

from copy import copy

import numpy as np
import pyarts3 as pyarts

POINTS = [(5e3, 10.0, 20.0), (42e3, -60.0, 170.0), (90e3, 75.0, -100.0)]

# %% Structs for pressure (hydrostatic), magnetic field (Schmidt-Legendre),
# wind (magnitude), and gravity (ellipsoid)

ws = pyarts.Workspace(False)
ws.surf_fieldPlanet(option="Earth")
ws.gravity_operatorCentralMass(mass=5.972168e24)
ws.atm_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.alt_grid = np.linspace(0, 100e3, 101)
ws.atm_fieldHydrostaticPressure(p0=1e5, fixed_specific_gas_constant=287.0)
for key, value in (("wind_u", 5.0), ("wind_v", -3.0), ("wind_w", 0.1)):
    wind = copy(ws.atm_field["t"])
    wind.data.data = np.full(wind.data.data.shape, value)
    ws.atm_field[key] = wind
ws.atm_fieldAbsoluteWindField()
ws.atm_fieldSchmidthFieldFromIGRF()

# %% The IGRF struct

igrf = pyarts.Workspace(False)
igrf.atm_fieldInit(toa=100e3, default_isotopologue="None")
igrf.atm_fieldIGRF()

KEYS = ["p", "wind_u", "wind_v", "wind_w", "mag_u", "mag_v", "mag_w"]

for filetype in ["ascii", "binary"]:
    ws.savexml("operators.xml", type=filetype)
    loaded = pyarts.Workspace.fromxml("operators.xml")

    igrf.atm_field.savexml("igrf.xml", type=filetype)
    loaded_igrf = pyarts.arts.AtmField.fromxml("igrf.xml")

    for point in POINTS:
        for key in KEYS:
            assert ws.atm_field[key].data(*point) == loaded.atm_field[key].data(*point), (
                filetype,
                key,
                point,
            )
        for key in ["mag_u", "mag_v", "mag_w"]:
            assert igrf.atm_field[key].data(*point) == loaded_igrf[key].data(*point), (
                filetype,
                key,
                point,
            )
        assert ws.gravity_operator(*point) == loaded.gravity_operator(*point), (
            filetype,
            point,
        )
