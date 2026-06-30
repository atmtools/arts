import numpy as np
import pyarts3 as pyarts

# x_ref: pinned strided values (stride 60000) of the native-grid
# spectral_propmat output for O3 at STP. Generated once from the
# native-grid case and frozen as the regression reference.
x_ref = [
    5.48830541e-29,
    7.01126943e-26,
    2.25566466e-25,
    3.67157516e-27,
    5.42679027e-26,
    7.56822053e-23,
    1.11608430e-21,
    3.90151892e-22,
]

# Tolerances: 1e-36 is the native-grid round-trip floor (effectively
# exact). The resampled grids go through np.linspace, which introduces
# interpolation error, so we loosen to 1e-18 for those cases.
cases = [
    (1.0, 1e-36),  # native grid; exact
    (0.2, 1e-18),  # downsampled; interpolation error dominates
    (1.2, 1e-18),  # upsampled;   same note as downsampled
]


def _check_grid(ws, atm, nd, f, stride, atol):
    x = ws.abs_xfit_data.spectral_propmat(f=f, atm=atm, spec="O3") / nd
    assert np.allclose(x[::stride, 0], x_ref, atol=atol)


def main():
    ws = pyarts.Workspace()

    atm = pyarts.arts.AtmPoint()
    atm.pressure = 101325.0
    atm.temperature = 273.15
    atm["O3"] = 1e-8
    nd = atm.number_density("O3-666")

    ws.abs_speciesSet(species=["O3-XFIT"])
    ws.ReadCatalogData()
    f = ws.abs_xfit_data["O3"].fitcoeffs[0].grids[0]

    for n, atol in cases:
        f_n = f if n == 1.0 else np.linspace(f[0], f[-1], int(f.size * n))
        _check_grid(ws, atm, nd, f_n, stride=int(60000 * n), atol=atol)


if __name__ == "__main__":
    main()
