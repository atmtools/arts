"""
Test calling of workspace methods.
"""
import os
import numpy as np
import pytest
import scipy as sp
import arts
from arts.workspace import Workspace
from arts.xml import load, save

class TestMethods:
    """
    Tests the calling of ARTS workspace methods.
    """
    def setup_method(self):
        """
        This ensures a new Workspace for every test.
        """
        self.ws  = Workspace(verbosity = 0)
        self.setup_workspace()

    def setup_workspace(self):
        self.ws.execute_controlfile("artscomponents/clearsky/TestClearSky.arts")

    def test_mixed_arguments(self):
        """
        Check that this raises a syntax error.
        """
        ws = self.ws
        with pytest.raises(SyntaxError):
            ws.yCalc(ws.yf, y_f = ws.y_f)

    def test_unexpected_argument(self):
        """
        Providing a named argument with a name that is not actually an argument
        of the WSM should raise an error.
        """
        ws = self.ws
        with pytest.raises(Exception):
            ws.yCalc(nonsense = ws.y_f)

    def test_override_output(self):
        """
        Test overriding of output parameters in the two possible
        ways.
        """
        ws = self.ws
        y_ref = np.copy(ws.y.value)

        ws.yf = np.zeros(y_ref.size)
        ws.yCalc(ws.yf)
        assert(np.allclose(ws.yf.value, y_ref))

        ws.yf = np.zeros(y_ref.size)
        ws.yCalc(y = ws.yf)
        assert(np.allclose(ws.yf.value, y_ref))

    def test_override_input(self):
        """
        Test overriding of input parameters in the two possible ways.

        atmgeom_checked WSV is set to zero so that both calculations
        should fail if input is not overridden.
        """
        ws = self.ws
        ws.atmgeom_checked = 0

        ws.rte_pos = np.array([600e3, 0, 0])
        ws.rte_pos2 = np.array([600e3, 0])
        ws.rte_los = np.array([180.0, 0])

        ws.ppathCalc(ws.ppath, ws.ppath_agenda, ws.ppath_lmax, ws.ppath_lraytrace, 1)
        ws.ppathCalc(atmgeom_checked = 1)

    def test_generic_input(self):
        """
        Test overriding of generic input in the two possible ways.
        """
        ws = self.ws
        species = (["H2O-SelfContStandardType, H2O-ForeignContStandardType, H2O",
                    "N2-SelfContStandardType",
                    "O3"] )

        ws.ArrayOfArrayOfSpeciesTagCreate("abs_species_2")
        ws.abs_speciesSet(ws.abs_species_2, ws.abs_xsec_agenda_checked,
                          ws.propmat_clearsky_agenda_checked, species)
        ws.ArrayOfArrayOfSpeciesTagCreate("abs_species_3")
        ws.abs_speciesSet(abs_species = ws.abs_species_3, species = species)
        assert(ws.abs_species_2.value == ws.abs_species_3.value)

    def test_generic_output(self):
        """
        Test overriding of generic input in the two possible ways.
        """
        ws = self.ws

        mat = np.ones((2, 2))
        save(mat, "matrix.xml")

        ws.sensor_los = np.zeros((2, 2))
        ws.ReadXML(ws.sensor_los, "matrix.xml")
        assert(np.allclose(mat, ws.sensor_los.value))

        ws.sensor_los = np.zeros((2, 2))
        ws.ReadXML(out = ws.sensor_los, filename = "matrix.xml")
        assert(mat, ws.sensor_los.value)

    def test_supergeneric_overload_resolution(self):
        """
        Test resolution of supergeneric methods.
        """
        self.ws.ArrayOfIndexCreate("array_of_index")
        self.ws.ArrayOfArrayOfIndexCreate("array_of_array_of_index")
        self.ws.array_of_index = [1, 2, 3]
        self.ws.Append(self.ws.array_of_array_of_index, self.ws.array_of_index)
        self.ws.Append(self.ws.array_of_array_of_index, self.ws.array_of_index)

    def test_supergeneric_overload_failure(self):
        """
        Test expected failure of supergeneric overload resolution.
        """
        with pytest.raises(Exception):
            self.ws.NumericCreate("numeric_wsv")
            self.ws.StringCreate("string_wsv")
            self.ws.Copy(self.ws.string_wsv, self.ws.numeric_wsv)

    def test_wsm_error(self):
        """
        Test error handling from ARTS WSMs.
        """
        with pytest.raises(Exception):
            self.ws.yCalc()
