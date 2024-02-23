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
    
    def test_empty(self):
        pass


if __name__ == "__main__":
    x = TestAgendas()
    x.test_empty()
