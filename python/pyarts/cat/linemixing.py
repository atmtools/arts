# -*- coding: utf-8 -*-

"""This contains methods to replace set line mixing
"""

import numpy as np

def init_ecs(ws, vmrs=np.array([0.21, 0.79]), specs=["O2", "N2"]):
    """ Initialize the error-corrected sudden molecular parameters
    
    Parameters:
        ws (Workspace): A pyarts workspace
        vmrs (ndarray): The volume mixing ration of "air"
        specs (str): A comma separated list of "air" species with resp. vmr
    
    Returns: None
    """
    ws.ecs_dataInit()
    ws.ecs_dataAddMakarov2020()
    ws.ecs_dataAddRodrigues1997()
    ws.ecs_dataAddTran2011()
    ws.ecs_dataSetMeanAir(vmrs=vmrs, specs=specs)

def adapt_lines(ws, lines, t_grid=np.linspace(150, 350),pressure=1e5,order=1,robust=1,rosenkranz_adaptation=0):
    """ Adapt the lines to use line mixing coefficients from on-the-fly bands
    
    Parameters:
        ws (Workspace): A pyarts workspace
        lines (pyarts workspace variable): An ArrayOfAbsorptionLines (modified in-place)
        t_grid (ndarray): A range of increasing temperatures in Kelvin
        pressure (float): A pressure in pascal
        order (int): The order of Rosenkranz adaptations
        robust (int): Whether or not an error is ignored
        rosenkranz_adaptation (int): Wheter or not pure perturbations [0] or eigenvalues are used in the adaptation [else]
        
    Returns:
        lines modified
    """
    ws.abs_linesAdaptOnTheFlyLineMixing(abs_lines=lines, t_grid=t_grid, pressure=pressure, order=order,
                                        robust=robust, rosenkranz_adaptation=rosenkranz_adaptation)
    return lines
