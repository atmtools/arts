import numpy as np
import pyarts

def main():
    solar_spec = pyarts.arts.GriddedField2()
    solar_spec.readxml("star/Sun/solar_spectrum_May_2004.xml")
    
    ws = pyarts.workspace.Workspace()
    ws.PlanetSet(option="Earth")
    
    ws.jacobianOff()
    
    ws.abs_speciesSet(species=['N2'])
    ws.abs_lines_per_speciesSetEmpty()
    ws.propmat_clearsky_agendaAuto()
    
    ws.stokes_dim = 1
    ws.f_grid = solar_spec.grids[0]
    ws.lat_true = [0.0]
    ws.lon_true = [0.0]
    ws.p_grid = np.logspace(5, -1, 2)
    ws.z_surface = [[300]]
    
    ws.AtmRawRead(basename="planets/Earth/Fascod/tropical/tropical")
    ws.AtmosphereSet1D()
    ws.AtmFieldsCalc()
    
    ws.surface_skin_t = ws.t_field.value[0, 0, 0]
    ws.surface_scalar_reflectivity = [0.3]
    
    ws.starsAddSingleFromGrid(
        star_spectrum_raw=solar_spec,
        latitude=0.0,
        longitude=0.0,
    )
    
    ws.sensorOff()
    ws.cloudboxSetFullAtm()
    ws.Touch(ws.scat_data)
    ws.pnd_fieldZero()
    
    
    @pyarts.workspace.arts_agenda(ws=ws, set_agenda=True)
    def gas_scattering_agenda(ws):
        ws.Ignore(ws.rtp_vmr)
        ws.gas_scattering_coefAirSimple()
        ws.gas_scattering_matRayleigh()
    
    
    ws.atmfields_checkedCalc()
    ws.lbl_checkedCalc()
    ws.scat_data_checkedCalc()
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    
    ws.DisortCalcIrradiance(nstreams=2)
    
    sim_irrad = (
        np.squeeze(ws.spectral_irradiance_field.value)[:, -1, 0] * -1
    )  # toa, downward
    data_irrad = flux_sun2flux_toa(ws, solar_spec.data[:, 0])  # stokes_dim 1
    
    assert np.allclose(
        sim_irrad, data_irrad
    ), "The simulated values do not match the solar ref spectra."


def flux_sun2flux_toa(ws, irrad):
    """scales the irradiance at Suns photosphere to TOA (100km)."""
    star_distance = ws.stars.value[0].distance
    star_radius = ws.stars.value[0].radius
    earth_radius = ws.refellipsoid.value[0]
    star_distance -= 100_000 + earth_radius

    factor = star_radius**2 / (star_radius**2 + star_distance**2)

    return irrad * factor

if __name__ == "__main__":
    main()
