"""
Test handling of workspace of the Python interface.
"""

import copy
import pytest
import pyarts3 as pyarts
import numpy as np
from pyarts3.workspace import Workspace, arts_agenda
from pyarts3.arts import Index, Vector


class TestWorkspace:
    def setup_method(self):
        pass

    def test_copy(self):
        ws = Workspace()

        ws.aaavar = Index(5)
        ws.aaavec = Vector([1, 2, 3])

        ws2 = copy.copy(ws)

        ws2.aaavar += 6
        ws2.aaavec += 10

        assert ws2.aaavar != ws.aaavar
        assert np.allclose(ws2.aaavec, ws.aaavec)

    def test_deepcopy(self):
        ws = Workspace()

        ws.aaavar = Index(5)
        ws.aaavec = Vector([1, 2, 3])

        ws2 = copy.deepcopy(ws)

        ws2.aaavar += 6
        ws2.aaavec += 10

        assert ws2.aaavar != ws.aaavar
        assert not np.allclose(ws2.aaavec, ws.aaavec)


if __name__ == "__main__":
    ta = TestWorkspace()
    ta.setup_method()
    ta.test_copy()
    ta.test_deepcopy()
