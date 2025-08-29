"""
Test calling of workspace methods.
"""
import numpy as np
import pytest
from tempfile import NamedTemporaryFile
import pyarts3 as pyarts
from pyarts3.utils.common import TempFileHandler
from pyarts3.workspace import Workspace


class TestMethods:
    """
    Tests the calling of ARTS workspace methods.
    """

    def setup_method(self):
        """
        This ensures a new Workspace for every test.
        """
        self.ws = Workspace()

    def test_unexpected_argument(self):
        """
        Providing a named argument with a name that is not actually an argument
        of the WSM should raise an error.
        """
        ws = self.ws
        with pytest.raises(Exception):
            ws.Print(nonsense=ws.y_f)

    def test_generic_input(self):
        """
        Test overriding of generic input in the two possible ways.
        """
        ws = self.ws
        species = [
            "H2O-SelfContStandardType, H2O-ForeignContStandardType, " "H2O",
            "N2-SelfContStandardType",
            "O3",
        ]
        
        ws.absorption_speciesSet(ws.absorption_species, species)
        ws.absorption_species_2 = pyarts.arts.ArrayOfArrayOfSpeciesTag()
        ws.absorption_speciesSet(ws.absorption_species_2, species)
        ws.absorption_species_3 = pyarts.arts.ArrayOfArrayOfSpeciesTag()
        ws.absorption_speciesSet(
            absorption_species=ws.absorption_species_3, species=species
        )
        assert ws.absorption_species == ws.absorption_species_3
        assert ws.absorption_species_2 == ws.absorption_species_3

    def test_generic_output(self):
        """
        Test overriding of generic input in the two possible ways.
        """
        ws = self.ws

        with TempFileHandler() as tempfile:
            ws.mat = pyarts.arts.Matrix()

            mat = np.ones((2, 2))
            ws.mat = np.ones((2, 2))
            ws.WriteXML("ascii", ws.mat, tempfile)

            ws.mat = np.zeros((2, 2))
            ws.ReadXML(ws.mat, tempfile)
            assert np.allclose(mat, ws.mat.value)

            ws.mat = np.zeros((2, 2))
            ws.ReadXML(output=ws.mat, filename=tempfile)
            assert np.allclose(mat, ws.mat)

    def test_supergeneric_overload_failure(self):
        """
        Test expected failure of supergeneric overload resolution.
        """
        self.ws.numeric_wsv = pyarts.arts.Numeric()
        self.ws.string_wsv = pyarts.arts.String()
        with pytest.raises(Exception):
            self.ws.Copy(self.ws.string_wsv, self.ws.numeric_wsv)

    def test_predefined_doc(self):
        m = "propagation_matrixAddPredefined"
        isots = pyarts.arts.globals.all_isotopologues()
        desc = pyarts.arts.globals.workspace_methods()[m].desc

        for x in isots:
            if x.predef:
                assert (
                    f"{x.name}:" in desc
                ), f"{x.name} is not documented properly (as '{x.name}: "
                f"SOME INFORMATION')\n\nPlease add this to {m}"


if __name__ == "__main__":
    x = TestMethods()
    x.setup_method()
    x.test_generic_output()
