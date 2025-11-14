"""
Test handling of agendas of the Python interface.
"""

import numpy as np
from pyarts3.workspace import Workspace, arts_agenda


def _input_only(freq_grid):
    assert np.allclose(freq_grid, [1, 2, 3])


def _return_only():
    freq_grid = [4, 5, 6]
    return freq_grid


def _input_output(freq_grid):
    assert np.allclose(freq_grid, [1, 2, 3])
    freq_grid = [4, 5, 6]
    return freq_grid


def _check_return_value(freq_grid):
    assert np.allclose(freq_grid, [4, 5, 6])


class TestAgendaPythonCalls:
    """
    Tests that test correct behaviour of calling Python
    functions inside an Agenda.
    """

    def init_spectral(self, ws):
        ws.spectral_radiance_observer_position = [1, 2, 3]
        ws.spectral_radiance_observer_line_of_sight = [1, 2]

    def test_python_function_return_only(self):
        ws = Workspace()

        @arts_agenda(ws=ws, fix=True)
        def ray_path_observer_agenda(ws):
            _return_only()
            _check_return_value()

        ws.ray_path_observer_agenda = ray_path_observer_agenda
        self.init_spectral(ws)

        ws.freq_grid = [1, 2, 3]
        ws.ray_path_observer_agendaExecute()
        assert np.allclose(ws.freq_grid, [1, 2, 3])

    def test_python_function_input_only(self):
        ws = Workspace()

        @arts_agenda(ws=ws, fix=True)
        def ray_path_observer_agenda(ws):
            _input_only()

        ws.ray_path_observer_agenda = ray_path_observer_agenda
        self.init_spectral(ws)

        ws.freq_grid = [1, 2, 3]
        ws.ray_path_observer_agendaExecute()
        assert np.allclose(ws.freq_grid, [1, 2, 3])

    def test_python_function_input_output(self):
        ws = Workspace()

        @arts_agenda(ws=ws, fix=True)
        def ray_path_observer_agenda(ws):
            _input_output()
            _check_return_value()

        ws.ray_path_observer_agenda = ray_path_observer_agenda
        self.init_spectral(ws)

        ws.freq_grid = [1, 2, 3]
        ws.ray_path_observer_agendaExecute()
        assert np.allclose(ws.freq_grid, [1, 2, 3])


if __name__ == "__main__":
    x = TestAgendaPythonCalls()
    x.test_python_function_input_only()
    x.test_python_function_input_output()
    x.test_python_function_return_only()
