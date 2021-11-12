# -*- coding: utf-8 -*-

import numpy as np

def select(ws, lines, qid, fmin=-np.infty, fmax=np.infty, safe=1):
    removed_lines = pyarts.classes.ArrayOfAbsorptionLines()
    removed_lines.set(lines)

    ws.abs_linesKeepBand(abs_lines=lines, qid=qid)
    ws.abs_linesRemoveEmptyBands(abs_lines=lines)
    ws.abs_linesRemoveLines(abs_lines=lines, lower_frequency=fmin, upper_frequency=fmax, safe=safe)

    assert len(lines) == 1, "Can only have one band"
    assert len(lines[0].lines), "Must have remaining lines"

    return lines
