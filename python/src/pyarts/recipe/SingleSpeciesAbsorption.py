import pyarts


class SingleSpeciesAbsorption:
    """Calculates absorption coefficients for a single absorbing species."""

    def __init__(
        self,
        species: str,
    ):
        """Initialization

        Parameters
        ----------
        species : str
            See absorption_speciesSet for details.
        """
        self.ws = pyarts.Workspace()
        self.ws.WignerInit()
        self.ws.absorption_speciesSet(species=[species])
        self.ws.ReadCatalogData()
        self.ws.propagation_matrix_agendaAuto()
        self.ws.ray_path_point = pyarts.arts.PropagationPathPoint()

    def __call__(
        self,
        frequency_grid: pyarts.arts.AscendingGrid,
        atmospheric_point: pyarts.arts.AtmPoint,
    ):
        """Call operator to return a propagation matrix

        Parameters
        ----------
        frequency_grid : ~pyarts.arts.AscendingGrid
            A list of frequency points.
        atmospheric_point : ~pyarts.arts.AtmPoint
            The state of the atmosphere at the point of interest

        Returns
        -------
        propagation_matrix : numpy.ndarray
            The propagation matrix at the frequency and point of interest
            Note that the first dimention is the size of the frequency
            grid and that the second dimension contains 7 variables, the
            first of which is unpolarized absorption.
        """

        self.ws.propagation_matrix_agendaExecute(
            frequency_grid=frequency_grid,
            atmospheric_point=atmospheric_point,
        )

        return 1.0 * self.ws.propagation_matrix[:, 0]


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    o2_spec = SingleSpeciesAbsorption(species="O2-66")
    f = np.linspace(1e12, 1000e12, 1000)

    atm = pyarts.arts.AtmPoint()
    atm.set_species("O2", 0.21)
    atm.set_species("H2O", 1e-3)
    atm.temperature = 273
    atm.pressure = 1e5
    plt.semilogy(f, o2_spec(f, atm))

    atm.temperature = 200
    atm.pressure = 1e4
    plt.semilogy(f, o2_spec(f, atm))
    atm.species("O2") * atm.isotopologue("O2-66")
