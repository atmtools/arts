# DISORT reference-test coverage

`DISOTESTAUX.f` is the canonical numerical reference.  A case is only marked
`ported` when CPP-DISORT is compared with its original flux and/or user-angle
radiance values; merely constructing the atmosphere or checking finite output
does not count.

| Problem | Cases | Current state | CPP limitation or required work |
|---|---|---|---|
| 1 | a–f | ported | All original user-angle radiances and fluxes pass. Exact-unity input uses the finite constant/linear conservative pair without changing the supplied albedo. |
| 2 | a–d | ported | All original user-angle radiances and fluxes pass, including exact-unity input. |
| 3 | a–b | ported | Arbitrary-angle IMS/TMS corrections are applied to the original user-angle radiances. |
| 4 | a–c | ported | All Haze-L fluxes and user-angle radiances are checked with arbitrary-angle IMS/TMS corrections. |
| 5 | a–b | ported | Cloud C.1 inputs and references are shared for reuse by scalar VDISORT tests. |
| 6 | a–h | ported | All active Fortran direct/diffuse/up/DFDT references are checked. The tests map Planck source polynomials, directional Kirchhoff emission, and Lambertian/Hapke reflection to the CPP API. The exactly non-scattering cases use `omega=1e-8` to regularize otherwise degenerate nonzero Fourier modes. |
| 7 | a–e | ported | All active Fortran direct/diffuse/up/DFDT references are checked, including the narrow-band small-signal case and the Hapke thermal boundary. |
| 8 | a–c | ported | All two-layer fluxes and four user-angle radiances are checked; inputs and references are shared for scalar VDISORT reuse. |
| 9 | a–c | ported | All multilayer fluxes and user-angle radiances are checked, including the fully general thermal/beam/Lambertian case. |
| 10 | a–b | ported | The formal `u_user()` solution at all four quadrature angles is compared with direct `u()` output; pointwise and bulk quadrature paths are also checked. |
| 11 | a–b | ported | The original homogeneous layer is compared with its equivalent three-layer subdivision for all radiances, fluxes, and `DFDT`. |
| 12 | a–b | ported physical equivalence | The complete one-layer solution is compared with the equivalent three-layer solution for corrected radiances, fluxes, and `DFDT`; the original absorption cutoff is intentionally not implemented. |
| 13 | a–d | ported physical output | The regular solves reproduce both shortcut albedo/transmission pairs. The `IBCND=1` shortcut itself is unavailable, but is not needed to test its physical results. |
| 14 | a–d | ported | Hapke, Cox-Munk, RPV, and Ross-Li raw BRDFs, their Fourier reconstruction, all user-angle radiances, fluxes, and `DFDT` are checked. CPP represents the exactly transparent layer by `tau=1e-12`, `omega=1e-8` to avoid a degenerate eigensystem. |
| 15 | a–d | ported | The Rayleigh/aerosol atmosphere uses all 600 Kokhanovsky moments for arbitrary-angle IMS/TMS corrections. Hapke, shadowed Cox-Munk, RPV, and Ross-Li cases are checked against every original radiance, flux, and `DFDT` value. |
| 16 | a | intentionally unsupported | CPP-DISORT is kept plane-parallel; pseudo-spherical direct-beam corrections are out of scope and are not planned. |
| 17 | a–b | ported | Delta-M accepts a general removed-peak Legendre profile, and the DISORT 4 Gaussian delta-M-plus helper reproduces all 1,080 original aerosol/cloud user-angle radiances. IMS/TMS is intentionally unavailable for this mode, matching DISORT 4.0.99. |

The older cDISORT suite has 48 subcases.  Its Problems 1–13 overlap the
Fortran configurations above.  cDISORT Problem 14b additionally calls the
separate `twostr()` solver; CPP-DISORT with `NQuad=2` is not an implementation
of that algorithm and must not be presented as an exact port.
