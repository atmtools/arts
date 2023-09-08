"""
Test handling of workspace of the Python interface.
"""
import copy
import pytest
import pyarts
from pyarts.workspace import Workspace, arts_agenda
from pyarts.arts import Index


class TestWorkspace:
    def setup_method(self):
        pass

    def test_copy(self):
        ws = Workspace()
        ws.aaavar = Index(5)

        ws2 = copy.copy(ws)
        ws2.aaavar = 6

        assert ws2.aaavar == ws.aaavar

    def test_deepcopy(self):
        ws = Workspace()
        ws.aaavar = Index(5)

        ws2 = copy.deepcopy(ws)
        ws2.aaavar = 6

        assert ws2.aaavar != ws.aaavar


if __name__ == "__main__":
    ta = TestWorkspace()
    ta.setup_method()
    ta.test_copy()
    ta.test_deepcopy()
