# -*- coding: utf-8 -*-

def read(ws, filename):
    """ Reads an online Hitran catalog of the latest ARTS specifications

    These specifications are:

    Field separator: [tab], Line endings: Windows (CR LF)

    .par_line
    qns'
    qns''

    Returns are reference to the generated abs_lines
    """
    ws.ReadHITRAN(filename=filename,
      globalquantumnumbers="v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 r Lambda ElectronState parity kronigParity S Sym alpha n tau",
      localquantumnumbers="J N Omega Ka Kc F",
      normalization_option="SFS",
      hitran_type="Online")
    ws.abs_linesRemoveUnusedLocalQuantumNumbers()

    return ws.abs_lines
