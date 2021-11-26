# -*- coding: utf-8 -*-

import numpy as np

def select(ws, lines, qid, fmin=-np.infty, fmax=np.infty, safe=1):
    """ Select a single band
    
    Parameters:
        ws (Workspace): A pyarts workspace
        lines (pyarts workspace variable): An ArrayOfAbsorptionLines (modified in-place)
        qid (pyarts workspace variable): An ArrayOfQuantumIdentity to select a band
        fmin (float): The minimum frequency
        fmax (float): The maximum frequency
        safe (int): Whether or not to disallow removing lines unsafely (e.g., with line mixing you need the whole band)
        
    Returns:
        lines modified
    """
    ws.abs_linesKeepBand(lines, qid)
    ws.abs_linesRemoveEmptyBands(lines)
    ws.abs_linesRemoveLines(lines, fmin, fmax, 0.0, safe)
    ws.abs_linesRemoveEmptyBands(lines)

    assert len(lines.value) == 1, "Can only have one band"
    assert len(lines.value[0].lines), "Must have remaining lines"

    return lines
