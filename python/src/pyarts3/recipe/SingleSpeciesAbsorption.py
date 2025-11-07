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
        frequency_grid : ~pyarts3.arts.AscendingGrid
            A list of frequency points.
        atmospheric_point : ~pyarts3.arts.AtmPoint
            The state of the atmosphere at the point of interest

        Returns
        -------
        numpy.ndarray : propagation_matrix
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
