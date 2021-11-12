# -*- coding: utf-8 -*-

import os


def set(ws, lines, directory):
    if "pyarts_zeeman_set_aoqi__" not in dir(ws):
        ws.create_variable("ArrayOfQuantumIdentifier", "pyarts_zeeman_set_aoqi__")

    if "pyarts_zeeman_set_v__" not in dir(ws):
        ws.create_variable("Vector", "pyarts_zeeman_set_v__")

    assert ws.pyarts_zeeman_set_aoqi__.group == "ArrayOfQuantumIdentifier"
    assert ws.pyarts_zeeman_set_v__.group == "Vector"

    for file in os.listdir(directory):
        if ".qid.xml" not in file: continue
        
        aoqi = os.path.join(directory, file)
        v = os.path.join(directory, file.replace(".qid.xml", ".g.xml"))

        assert os.path.exists(aoqi), f"Cannot find file {aoqi}"  # How? I just asked for it!
        assert os.path.exists(v), f"Cannot find file {v}"  # We don't yet know if it exists!

        ws.ReadXML(ws.pyarts_zeeman_set_aoqi__, aoqi)
        ws.ReadXML(ws.pyarts_zeeman_set_v__, os.path.join(directory, v))
        ws.abs_linesSetZeemanCoefficients(lines, ws.pyarts_zeeman_set_aoqi__, ws.pyarts_zeeman_set_v__)
    
    return lines
