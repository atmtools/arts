# -*- coding: utf-8 -*-

"""This contains methods to replace set line mixing
"""

import numpy as np

def init_ecs(ws, vmrs=np.array([0.21, 0.79]), specs="O2, N2"):
    ws.ecs_dataInit()
    ws.ecs_dataAddMakarov2020()
    ws.ecs_dataAddRodrigues1997()
    ws.ecs_dataAddTran2011()
    ws.ecs_dataSetMeanAir(vmrs=vmrs, specs=specs)

def adapt_lines(ws, lines, t_grid=np.linspace(150, 350),pressure=1e5,order=1,robust=1,rosenkranz_adaptation=0):
    ws.abs_linesAdaptOnTheFlyLineMixing(abs_lines=lines, t_grid=t_grid, pressure=pressure, order=order,
                                        robust=robust, rosenkranz_adaptation=rosenkranz_adaptation)
    return lines
