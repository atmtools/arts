# -*- coding: utf-8 -*-

import os


def set(ws, lines, directory, spec):
    """ Set Zeeman coefficients
    
    Looks in directory for a file {spec}.aoqi.xml and {spec}.g.xml to find
    replacement data
    
    Parameters:
        ws (Workspace): A pyarts workspace
        lines (pyarts workspace variable): An ArrayOfAbsorptionLines (modified in-place)
        directory (str): A directory with Zeeman data
        spec (str): An Arts species
        
    Returns:
        lines modified
    """
    
    if "pyarts_zeeman_set_aoqi__" not in dir(ws):
        ws.create_variable("ArrayOfQuantumIdentifier", "pyarts_zeeman_set_aoqi__")

    if "pyarts_zeeman_set_v__" not in dir(ws):
        ws.create_variable("Vector", "pyarts_zeeman_set_v__")

    assert ws.pyarts_zeeman_set_aoqi__.group == "ArrayOfQuantumIdentifier"
    assert ws.pyarts_zeeman_set_v__.group == "Vector"
    
    aoqi = os.path.join(directory, f"{spec}.qid.xml")
    v = os.path.join(directory, f"{spec}.g.xml")

    assert os.path.exists(aoqi), f"Cannot find file {aoqi}"
    assert os.path.exists(v), f"Cannot find file {v}"

    ws.ReadXML(ws.pyarts_zeeman_set_aoqi__, aoqi)
    ws.ReadXML(ws.pyarts_zeeman_set_v__, os.path.join(directory, v))
    ws.abs_linesSetZeemanCoefficients(lines, ws.pyarts_zeeman_set_aoqi__, ws.pyarts_zeeman_set_v__)
    
    return lines
