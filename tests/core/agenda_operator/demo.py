import pyarts3 as pyarts
import numpy as np


def interface_function(freq, wind_jac, jac, species, path_point, atm):
    """ The interface function takes the same input as the agenda, in the
    same order.  The output must also match the agenda output.

    To see the input and outputs of an Agenda, write ws.agenda? in an ipython
    terminal or check in the pyarts documentation - especially the
    call-operator of agendaOperator classes describe in full the interface

    Anyways, the agenda spectral_propmat_agenda takes:
        a frequency grid
        a wind-shift-jacobian
        the jacobian target lists
        a selection of absorption species
        a ray path point
        the atmospheric state

    it outputs:
        propmat
        source-vec
        propmat-jacobian
        source-vec-jacobian

    The below demonstrates an incredibly easy way to set up these calculations.

    It is not intent for anything but documenting the idea, you need to fill
    this up yourself!
    """
    ws = pyarts.Workspace()
    ws.spectral_propmatInit(freq_grid=freq, jac_targets=jac)
    ws.spectral_propmat[:, 0] = 1
    return (
        ws.spectral_propmat,
        ws.spectral_nlte_srcvec,
        ws.spectral_propmat_jac,
        ws.spectral_nlte_srcvec_jac,
    )


ws = pyarts.Workspace()
ws.spectral_propmat_agendaSetOperator(f=interface_function)

ws.freq_grid = [1, 2]
ws.freq_wind_shift_jac
ws.ray_point
ws.atm_point

ws.spectral_propmat_agendaExecute()

assert np.allclose(
    ws.spectral_propmat, [[1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0]]
)
