#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 16:19:05 2021

@author: u237023
"""

import numpy as np
import pyarts

ws = pyarts.workspace.Workspace()

ws.isotopologue_ratiosInitFromBuiltin()
ws.Wigner6Init(ws.wigner_initialized)

ws.Touch(ws.jacobian_quantities)
ws.Touch(ws.rtp_nlte)

ws.abs_speciesSet(species=["O2-66"])
ws.ReadXML(ws.abs_lines, "lines/O2-66.xml")

qn = pyarts.classes.QuantumIdentifier("O2-66 S 1 1 Lambda 0 0 v1 0 0 ElecStateLabel X X")

ws.abs_linesKeepBand(qid=qn)
ws.abs_linesRemoveEmptyBands()
ws.abs_linesRemoveLines(lower_frequency=30e9, upper_frequency=120e9, safe=0)

assert len(ws.abs_lines.value) == 1

ws.abs_lines_per_speciesCreateFromLines()

ws.propmat_clearsky_agenda_checked = 1
ws.lbl_checked = 1
ws.stokes_dim = 1
ws.nlte_do = 0

ws.f_grid = np.linspace(40e9, 130e9, 1001)
ws.rtp_pressure = 1e5
ws.rtp_temperature = 273
ws.rtp_vmr = np.array([0.21])

ws.propmat_clearskyInit()
ws.propmat_clearskyAddLines()

old = 1.0 * ws.propmat_clearsky.value.data.value.flatten()

ws.ecs_dataInit()
ws.ecs_dataAddRodrigues1997()
ws.ecs_dataAddTran2011()
ws.ecs_dataAddMakarov2020()
ws.ecs_dataSetMeanAir(vmrs=[1], specs="N2")  # Since the rest is for O2

ws.abs_lines_per_speciesSetPopulation(option="ByMakarovFullRelmat")
ws.abs_lines_per_speciesAdaptOnTheFlyLineMixing(t_grid = np.linspace(150, 350, 51), pressure=1e5, order=1)

new = 1.0 * ws.propmat_clearsky.value.data.value.flatten()

assert np.isclose(old, new).all()
