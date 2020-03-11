# Included by seviri_fast.arts and seviri_reference.arts

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.abs_speciesSet(
    species=[
        "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252",
        "O3",
        "CO2, CO2-CKDMT252",
        "N2O",
        "CO",
        "CH4",
        "O2, O2-CIAfunCKDMT100",
        "N2, N2-CIAfunCKDMT252, N2-CIArotCKDMT252",
    ]
)
# abs_lineshape_per_tgDefine( abs_lineshape,
#                             abs_species,
#                             ["Voigt_Kuntz6", "Voigt_Kuntz6", "Voigt_Kuntz6", "Voigt_Kuntz6",
# 			     "Voigt_Kuntz6", "Voigt_Kuntz6", "Voigt_Kuntz6", "Voigt_Kuntz6"],
#                             ["VVH", "VVH", "VVH", "VVH", "VVH", "VVH", "VVH", "VVH"],
#                             [750e9, 750e9, 750e9, 750e9, 750e9, 750e9, 750e9, 750e9] )
ws.abs_lineshapeDefine(ws.abs_lineshape, "Voigt_Kuntz6", "VVH", 750000000000.0)
