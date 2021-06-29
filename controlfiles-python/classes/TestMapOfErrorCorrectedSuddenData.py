#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 16:47:44 2021

@author: larsson
"""

import pyarts
import numpy as np
import matplotlib.pyplot as plt


arts = pyarts.workspace.Workspace()
arts.isotopologue_ratiosInitFromBuiltin()

arts.ecs_dataInit()
arts.ecs_dataAddRodrigues1997()
arts.ecs_dataAddMakarov2020()
arts.ecs_data.value.savexml("tmp.xml")

ecs_data2 = pyarts.classes.MapOfErrorCorrectedSuddenData()
ecs_data2.readxml("tmp.xml")

# Cannot compare directly but this gets close enough to ensure reading works
assert np.isclose(float(ecs_data2["CO2-626 ALL"]["O2"].mass),
                  float(arts.ecs_data.value["CO2-626 ALL"]["O2"].mass))
