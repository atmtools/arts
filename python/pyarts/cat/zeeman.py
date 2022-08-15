# -*- coding: utf-8 -*-

import os
import pyarts


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
    
    aoqi = pyarts.arts.ArrayOfQuantumIdentifier()
    vec = pyarts.arts.Vector()
    
    aoqi_file = os.path.join(directory, f"{spec}.qid.xml")
    vec_file = os.path.join(directory, f"{spec}.g.xml")

    assert os.path.exists(aoqi_file), f"Cannot find file {aoqi_file}"
    assert os.path.exists(vec_file), f"Cannot find file {vec_file}"

    aoqi.readxml(aoqi_file)
    vec.readxml(vec_file)
    ws.abs_linesZeemanCoefficients(lines, aoqi, vec)
    
    return lines
