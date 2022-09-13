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
        assert ws.aaavar.value == ws.aaavar.value
        
        ws2 = copy.copy(ws)
        ws2.aaavar = 6
        assert ws2.aaavar.value == ws2.aaavar.value
        
        assert ws2.aaavar.value == ws.aaavar.value
    
    def test_deepcopy(self):
        ws = Workspace()
        ws.aaavar = Index(5)
        assert ws.aaavar.value == ws.aaavar.value
        
        ws2 = copy.deepcopy(ws)
        ws2.aaavar = 6
        assert ws2.aaavar.value == ws2.aaavar.value
        
        assert ws2.aaavar.value != ws.aaavar.value
    
    def test_copy_agenda(self):
        ws = Workspace()
        ws.aaavar = Index(5)
        
        @arts_agenda(ws=ws, set_agenda=True)
        def test_agenda(ws):
            ws.Print(ws.aaavar, 0)
        ws.AgendaExecute(ws.test_agenda)
        
        ws2 = copy.deepcopy(ws)
        ws2.aaavar = 4
        ws2.AgendaExecute(ws2.test_agenda)
        

if __name__ == "__main__":
    ta = TestWorkspace()
    ta.setup_method()
    ta.test_copy()
    ta.test_deepcopy()
    ta.test_copy_agenda()
    