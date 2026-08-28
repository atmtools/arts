# DISORT reference-test coverage

`DISOTESTAUX.f` is the canonical numerical reference.  A case is only marked
`ported` when CPP-DISORT is compared with its original flux and/or user-angle
radiance values; merely constructing the atmosphere or checking finite output
does not count.

| Problem | Cases | Current state | CPP limitation or required work |
|---|---|---|---|
| 1 | a–f | ported | All original user-angle radiances and fluxes pass, including exact-unity input through the internal conservative dither. |
| 2 | a–d | ported | All original user-angle radiances and fluxes pass, including exact-unity input. |
| 3 | a–b | ported | Arbitrary-angle IMS/TMS corrections are applied to the original user-angle radiances. |
| 4 | a–c | ported | All Haze-L fluxes and user-angle radiances are checked with arbitrary-angle IMS/TMS corrections. |
| 5 | a–b | ported | Cloud C.1 inputs and references are shared for reuse by scalar VDISORT tests. |
| 6 | a–h | smoke only; references extracted | Active Fortran direct/diffuse/up/DFDT values are now shared; boundary and thermal-source inputs still require normalization mapping to the CPP API. |
| 7 | a–e | smoke only; references extracted | Active Fortran direct/diffuse/up/DFDT values are now shared; numerical assertions require the thermal/boundary mapping. |
| 8 | a–c | ported | All two-layer fluxes and four user-angle radiances are checked; inputs and references are shared for scalar VDISORT reuse. |
| 9 | a–c | ported | All multilayer fluxes and user-angle radiances are checked, including the fully general thermal/beam/Lambertian case. |
| 10 | a–b | ported | The formal `u_user()` solution at all four quadrature angles is compared with direct `u()` output; pointwise and bulk quadrature paths are also checked. |
| 11 | a–b | ported | The original homogeneous layer is compared with its equivalent three-layer subdivision for all radiances, fluxes, and `DFDT`. |
| 12 | a–b | ported physical equivalence | The complete one-layer solution is compared with the equivalent three-layer solution for corrected radiances, fluxes, and `DFDT`; the original absorption cutoff is intentionally not implemented. |
| 13 | a–d | references extracted | CPP lacks the `IBCND=1` shortcut; its two albedo/transmission pairs are shared and can be tested through the regular method. |
| 14 | a–d | pending | Fortran 4.0.99 BRDF cases are portable through CPP BRDF callbacks. |
| 15 | a–d | pending | Multilayer BRDF cases are portable through CPP BRDF callbacks. |
| 16 | a | blocked by solver feature | CPP has no pseudo-spherical direct-beam correction. |
| 17 | a–b | blocked by solver feature | CPP has delta-M scaling, but not DISORT 4.0.99 delta-M-plus. |

The older cDISORT suite has 48 subcases.  Its Problems 1–13 overlap the
Fortran configurations above.  cDISORT Problem 14b additionally calls the
separate `twostr()` solver; CPP-DISORT with `NQuad=2` is not an implementation
of that algorithm and must not be presented as an exact port.
