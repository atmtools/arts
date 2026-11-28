import numpy as np
import pyarts3 as pyarts

A = pyarts.arts

# A callback gets the point and the spectral point of whichever medium it is called in, and deals with the
# media itself
seen = []


def bulk(point, spectral_point, f, n):
    seen.append((type(point).__name__, type(spectral_point).__name__))
    if isinstance(point, A.SubsurfacePoint):
        raise RuntimeError("this callback has no subsurface optics")
    extinction = np.asarray(spectral_point["extinction"])
    return A.ScatteringTroSpectralVector(
        phase_matrix=A.CompactPlanarMuelmatMatrix(np.zeros((len(f), n + 1, 6))),
        extinction_matrix=A.CompactPlanarPropmatVector(np.column_stack([extinction, np.zeros((len(f), 2))])),
        absorption_vector=A.CompactPlanarStokvecVector(np.column_stack([extinction, np.zeros(len(f))])),
    )


atm_spectral = A.AtmSpectralPoint()
atm_spectral["extinction"] = [1e-3, 2e-3]

species = A.CallbackScatteringSpecies(A.ScatteringGeneralSpectralTRO(bulk))
out = species.tro_spectral(A.AtmPoint(), atm_spectral, [1e9, 2e9], 3)
assert np.asarray(out.phase_matrix).shape == (2, 4, 6)
assert np.array_equal(np.asarray(out.extinction_matrix)[:, 0], [1e-3, 2e-3])

try:
    species.tro_spectral(A.SubsurfacePoint(), A.SubsurfaceSpectralPoint(), [1e9], 3)
    raise AssertionError("the callback's own error must reach the caller")
except AssertionError:
    raise
except Exception as e:
    assert "no subsurface optics" in str(e), str(e)

assert seen == [("AtmPoint", "AtmSpectralPoint"), ("SubsurfacePoint", "SubsurfaceSpectralPoint")], seen

# A point and a spectral point of different media are no call
try:
    species.tro_spectral(A.AtmPoint(), A.SubsurfaceSpectralPoint(), [1e9], 3)
    raise AssertionError("an atmospheric point with a subsurface spectral point must fail")
except TypeError:
    pass

# What the callback gives must have the shape asked for
wrong = A.CallbackScatteringSpecies(A.ScatteringGeneralSpectralTRO(lambda p, s, f, n: bulk(p, s, f, n + 1)))
one = A.AtmSpectralPoint()
one["extinction"] = [1e-3]
try:
    wrong.tro_spectral(A.AtmPoint(), one, [1e9], 3)
    raise AssertionError("a series of the wrong degree must fail")
except AssertionError:
    raise
except Exception as e:
    assert "Legendre series" in str(e), str(e)

# A callback is a scattering species
array = A.ArrayOfScatteringSpecies([A.ScatteringGeneralSpectralTRO(bulk)])
assert np.asarray(array.tro_spectral(A.AtmPoint(), one, [1e9], 2).phase_matrix).shape == (1, 3, 6)
