import pyarts
import numpy as np
import os


def create_local(gfn, lfn):
    """ Simple helper method to regenerate this data in case necessary """
    assert gfn != lfn, "Please use separate filenames"

    with open(gfn, "r") as a:
        line = a.readline()
        single_lines = {}
        while line != "":
            spec = line[:3]
            if spec not in single_lines:
                single_lines[spec] = line
            line = a.readline()

    with open(lfn, "w") as a:
        for isot in single_lines:
            a.write(single_lines[isot])


DIR = "lines"
local_fn = "single_line.par"

ws = pyarts.Workspace()
ws.absorption_bandsReadHITRAN(file="/Users/richard/Work/hitran.par.txt")
ws.absorption_bandsSaveSplit(dir=DIR)

# create_local("a-full-hitran-to-extract-from.par", local_fn)

# ws = pyarts.Workspace()

# ws.absorption_bandsReadHITRAN(file=local_fn)
# tmp = pyarts.arts.AbsorptionBands(ws.absorption_bands)
# ws.absorption_bandsSaveSplit(dir=DIR)
# ws.absorption_bands = {}
# ws.absorption_bandsReadSplit(dir=DIR)

# simple_attrs = ["a", "f0", "e0", "gu", "gl"]

# for key in tmp:
#     assert key in ws.absorption_bands, "Missing isotopologue key after read"

#     assert len(ws.absorption_bands[key].lines) == len(tmp[key].lines)
#     assert len(tmp[key].lines) == 1

#     l1 = ws.absorption_bands[key].lines[0]
#     l2 = tmp[key].lines[0]

#     # Check pure numeric values are reproduced
#     for attr in simple_attrs:
#         assert np.allclose(getattr(l1, attr), getattr(l2, attr)), (
#             f"Attribute '{attr}' not close "
#             + f"l1.{attr}={getattr(l1, attr)} vs l2.{attr}={getattr(l2, attr)}"
#         )

#     # Check Zeeman effect is reproduced
#     assert np.allclose(l1.z.gl, l2.z.gl), "Bad Zeeman GL"
#     assert np.allclose(l1.z.gu, l2.z.gu), "Bad Zeeman GU"
#     assert l1.z.on == l2.z.on

#     # Check that the line shape is correct
#     assert np.allclose(l1.ls.one_by_one, l2.ls.one_by_one)
#     assert np.allclose(l1.ls.T0, l2.ls.T0)
#     assert len(l1.ls.single_models) == len(l2.ls.single_models)
#     for i in range(len(l1.ls.single_models)):
#         lsm1 = l1.ls.single_models[i]
#         lsm2 = l2.ls.single_models[i]
#         assert lsm1.species == lsm2.species
#         assert len(lsm1.data) == len(lsm2.data)

#         for j in range(len(lsm1.data)):
#             lsmd1 = lsm1.data[j]
#             lsmd2 = lsm2.data[j]
#             assert lsmd1[0] == lsmd2[0]
#             assert lsmd1[1].type == lsmd2[1].type
#             assert np.allclose(lsmd1[1].data, lsmd2[1].data)

#     # Check quantum numbers
#     qn1 = l1.qn
#     qn2 = l2.qn
#     assert len(qn1) == len(qn2)
#     for q1, q2 in zip(qn1, qn2):
#         assert q1.type == q2.type
#         assert q1.low == q2.low
#         assert q1.upp == q2.upp

# # CLEAN-UP on success
# for fn in os.listdir(DIR):
#     path = os.path.join(DIR, fn)
#     os.remove(path)
# os.rmdir(DIR)
