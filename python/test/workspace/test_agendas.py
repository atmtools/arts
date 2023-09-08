"""
Test handling of agendas of the Python interface.
"""
import os
import numpy as np
import pytest
import scipy as sp
import pyarts
from pyarts.workspace import Workspace, arts_agenda
from pyarts.xml import load, save


class TestAgendas:
    """
    Tests the calling of ARTS workspace methods.
    """
    
    def test_planet_set(self):
        options = pyarts.arts.options.planetDefaultOptions.get_options_as_strings()

        ws = pyarts.workspace.Workspace()
        for opt in options:
            assert opt in ws.PlanetSet.__doc__, f"The {opt}-option is not documented correctly"
            ws.PlanetSet(option=opt)


if __name__ == "__main__":
    pass

