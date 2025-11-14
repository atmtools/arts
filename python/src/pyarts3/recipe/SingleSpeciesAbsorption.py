import pyarts3 as pyarts


class SingleSpeciesAbsorption:
    """Calculates absorption coefficients for a single absorbing species."""

    def __init__(
        self,
        species: str,
        cutoff: float = None,
    ):
        """Initialization

        Parameters
        ----------
        species : str
            See abs_speciesSet for details.
        cutoff : float
            The cutoff value for the absorption bands. Defaults to None for no cutoff.
        """
        self.ws = pyarts.Workspace()
        self.ws.WignerInit()
        self.ws.abs_speciesSet(species=[species])
        self.ws.ReadCatalogData()
        if cutoff is not None:
            for band in self.ws.abs_bands:
                self.ws.abs_bands[band].cutoff = "ByLine"
                self.ws.abs_bands[band].cutoff_value = cutoff
        self.ws.spectral_propmat_agendaAuto()
        self.ws.ray_path_point = pyarts.arts.PropagationPathPoint()

    def __call__(
        self,
        freq_grid: pyarts.arts.AscendingGrid,
        atm_point: pyarts.arts.AtmPoint,
    ):
        """Call operator to return a propagation matrix

        Parameters
        ----------
        freq_grid : ~pyarts3.arts.AscendingGrid
            A list of frequency points.
        atm_point : ~pyarts3.arts.AtmPoint
            The state of the atmosphere at the point of interest

        Returns
        -------
        numpy.ndarray : spectral_propmat
            The propagation matrix at the frequency and point of interest
            Note that the first dimention is the size of the frequency
            grid and that the second dimension contains 7 variables, the
            first of which is unpolarized absorption.
        """

        self.ws.spectral_propmat_agendaExecute(
            freq_grid=freq_grid,
            atm_point=atm_point,
        )

        return 1.0 * self.ws.spectral_propmat[:, 0]
