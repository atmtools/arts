import pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

# %% Sensor

ws.frequency_grid = np.linspace(-50e6, 50e6, 11) + 118750348044.712

# %% Species and line absorption

ws.absorption_speciesSet(species=["O2-66-118e9-119e9"])
ws.ReadCatalogData()
ws.absorption_bandsSetZeeman(isot="O2-66", fmin=118e9, fmax=119e9)
ws.WignerInit()

# %% Use the automatic agenda setter for propagation matrix calculations
ws.propagation_matrix_agendaAuto()

# %% Grids and planet

ws.surface_fieldSetPlanetEllipsoid(option="Earth")
ws.surface_field[pyarts.arts.options.SurfaceKey("t")] = 295.0
ws.atmospheric_fieldRead(
    toa=100e3, basename="planets/Earth/afgl/tropical/", missing_is_zero=1
)
ws.atmospheric_fieldIGRF(time="2000-03-11 14:39:37")

# %% Agenda settings

ws.spectral_radiance_space_agendaSet(option="UniformCosmicBackground")
ws.spectral_radiance_surface_agendaSet(option="Blackbody")
ws.spectral_radiance_observer_agendaSet(option="Emission")
ws.propagation_path_observer_agendaSet(option="Geometric")

# %% Main computations

sensor_pos = [[200e3, 0, 0], [300e3, 10, 0], [400e3, 20, 0], [500e3, 30, 0]]
sensor_los = [[180, 0], [179, 0], [178, 0], [177, 0]]

sensor_radiance = pyarts.arts.StokvecVector()
sensor_radiance_jacobian = pyarts.arts.StokvecMatrix()

ws.sensor_radianceFromObservers(
    sensor_radiance, sensor_radiance_jacobian, pos=sensor_pos, los=sensor_los
)

assert np.allclose(
    sensor_radiance.flatten(),
    np.array(
        [
            9.74014749e-16,
            -5.73266410e-22,
            3.35971333e-22,
            5.79543405e-20,
            9.87508131e-16,
            -8.85468548e-22,
            5.20680402e-22,
            7.15875974e-20,
            1.00474663e-15,
            -1.56097031e-21,
            9.22177956e-22,
            9.47601221e-20,
            1.02897295e-15,
            -3.46628874e-21,
            2.06369175e-21,
            1.41704789e-19,
            1.06996131e-15,
            -1.25441766e-20,
            7.61113565e-21,
            2.71627943e-19,
            8.92464529e-16,
            -2.74835774e-17,
            9.40358944e-18,
            1.64740750e-23,
            1.07031741e-15,
            -1.25526729e-20,
            7.61635272e-21,
            -2.71796673e-19,
            1.02965783e-15,
            -3.47094222e-21,
            2.06647677e-21,
            -1.41877124e-19,
            1.00574988e-15,
            -1.56411752e-21,
            9.24042794e-22,
            -9.49324644e-20,
            9.88822991e-16,
            -8.87855730e-22,
            5.22086678e-22,
            -7.17608757e-20,
            9.75636029e-16,
            -5.75202732e-22,
            3.37107304e-22,
            -5.81289199e-20,
        ]
    ),
)
