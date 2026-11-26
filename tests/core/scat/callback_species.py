import numpy as np
import pyarts3 as pyarts

A = pyarts.arts

# A callback gets the point of whichever medium it is called in, and deals with the media itself
seen = []


def bulk(point, f, n):
    seen.append(type(point).__name__)
    if isinstance(point, A.SubsurfacePoint):
        raise RuntimeError("this callback has no subsurface optics")
    return A.ScatteringTroSpectralVector(
        phase_matrix=A.CompactPlanarMuelmatMatrix(np.zeros((len(f), n + 1, 6))),
        extinction_matrix=A.CompactPlanarPropmatVector(np.zeros((len(f), 3))),
        absorption_vector=A.CompactPlanarStokvecVector(np.zeros((len(f), 2))),
    )


species = A.CallbackScatteringSpecies(A.ScatteringGeneralSpectralTRO(bulk))
out = species.tro_spectral(A.AtmPoint(), [1e9, 2e9], 3)
assert np.asarray(out.phase_matrix).shape == (2, 4, 6)

try:
    species.tro_spectral(A.SubsurfacePoint(), [1e9], 3)
    raise AssertionError("the callback's own error must reach the caller")
except AssertionError:
    raise
except Exception as e:
    assert "no subsurface optics" in str(e), str(e)

assert seen == ["AtmPoint", "SubsurfacePoint"], seen

# What the callback gives must have the shape asked for
wrong = A.CallbackScatteringSpecies(A.ScatteringGeneralSpectralTRO(lambda p, f, n: bulk(p, f, n + 1)))
try:
    wrong.tro_spectral(A.AtmPoint(), [1e9], 3)
    raise AssertionError("a series of the wrong degree must fail")
except AssertionError:
    raise
except Exception as e:
    assert "Legendre series" in str(e), str(e)

# A callback is a scattering species
array = A.ArrayOfScatteringSpecies([A.ScatteringGeneralSpectralTRO(bulk)])
assert np.asarray(array.tro_spectral(A.AtmPoint(), [1e9], 2).phase_matrix).shape == (1, 3, 6)
