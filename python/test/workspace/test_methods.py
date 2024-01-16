"""
Test calling of workspace methods.
"""
import numpy as np
import pytest
from tempfile import NamedTemporaryFile
import pyarts
from pyarts.workspace import Workspace


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

        ws.abs_speciesSet(
            ws.abs_species, ws.propmat_clearsky_agenda_checked, species
        )
        ws.abs_species_2 = pyarts.arts.ArrayOfArrayOfSpeciesTag()
        ws.abs_speciesSet(
            ws.abs_species_2, ws.propmat_clearsky_agenda_checked, species
        )
        ws.abs_species_3 = pyarts.arts.ArrayOfArrayOfSpeciesTag()
        ws.abs_speciesSet(abs_species=ws.abs_species_3, species=species)
        assert ws.abs_species == ws.abs_species_3
        assert ws.abs_species_2 == ws.abs_species_3

    def test_generic_output(self):
        """
        Test overriding of generic input in the two possible ways.
        """
        ws = self.ws

        tempfile = NamedTemporaryFile()

        mat = np.ones((2, 2))
        ws.avk = np.ones((2, 2))
        ws.WriteXML("ascii", ws.avk, tempfile.name)

        ws.avk = np.zeros((2, 2))
        ws.ReadXML(ws.avk, tempfile.name)
        assert np.allclose(mat, ws.avk.value)

        ws.avk = np.zeros((2, 2))
        ws.ReadXML(output=ws.avk, filename=tempfile.name)
        assert np.allclose(mat, ws.avk.value)

    def test_supergeneric_overload_failure(self):
        """
        Test expected failure of supergeneric overload resolution.
        """
        self.ws.numeric_wsv = pyarts.arts.Numeric()
        self.ws.string_wsv = pyarts.arts.String()
        with pytest.raises(Exception):
            self.ws.Copy(self.ws.string_wsv, self.ws.numeric_wsv)

    def test_predefined_doc(self):
        m = "propmat_clearskyAddPredefined"
        isots = pyarts.arts.globals.get_isotopologues()
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
    x.test_predefined_doc()
