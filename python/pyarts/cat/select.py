# -*- coding: utf-8 -*-

import numpy as np

def select(ws, lines, qid, fmin=-np.infty, fmax=np.infty, safe=1):
    ws.abs_linesKeepBand(lines, qid)
    ws.abs_linesRemoveEmptyBands(lines)
    ws.abs_linesRemoveLines(lines, fmin, fmax, 0.0, safe)
    ws.abs_linesRemoveEmptyBands(lines)

    assert len(lines.value) == 1, "Can only have one band"
    assert len(lines.value[0].lines), "Must have remaining lines"

    return lines
