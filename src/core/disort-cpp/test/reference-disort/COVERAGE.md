# DISORT reference-test coverage

`DISOTESTAUX.f` is the canonical numerical reference.  A case is only marked
`ported` when CPP-DISORT is compared with its original flux and/or user-angle
radiance values; merely constructing the atmosphere or checking finite output
does not count.

| Problem | Cases | Current state | CPP limitation or required work |
|---|---|---|---|
| 1 | a–f | ported | All original user-angle radiances and fluxes pass. Exact conservative scattering in 1b/1e is evaluated as `1-epsilon`, matching the existing CPP limitation. |
| 2 | a–d | ported | All original user-angle radiances and fluxes pass. Conservative cases use `1-epsilon`. |
| 3 | a–b | fluxes portable; user radiances blocked | Delta-scaled original radiances require IMS/TMS corrections evaluated at arbitrary user angles. |
| 4 | a–c | fluxes portable; user radiances blocked | Higher phase moments trigger delta scaling; arbitrary-angle IMS/TMS correction is incomplete. |
| 5 | a–b | fluxes portable; user radiances blocked | Delta-M intensity corrections at user angles remain incomplete. |
| 6 | a–h | smoke only | Boundary and thermal-source inputs require normalization mapping to the CPP API. |
| 7 | a–e | smoke only | Directly portable after thermal/boundary mapping. |
| 8 | a–c | pending | Directly portable. |
| 9 | a–c | 9c smoke only | Directly portable; general thermal, multilayer, beam, and surface case. |
| 10 | a–b | equivalence smoke only | CPP has formal user angles and quadrature output, so both paths are portable. |
| 11 | a–b | Pythonic references only | Original single-layer/multilayer values remain to be added. |
| 12 | a–b | smoke only | CPP lacks the absorption-optical-depth shortcut switch; physical output can be tested, shortcut execution cannot. |
| 13 | a–d | smoke only | CPP lacks the `IBCND=1` albedo/transmission shortcut; regular-method outputs can be tested. |
| 14 | a–d | pending | Fortran 4.0.99 BRDF cases are portable through CPP BRDF callbacks. |
| 15 | a–d | pending | Multilayer BRDF cases are portable through CPP BRDF callbacks. |
| 16 | a | blocked by solver feature | CPP has no pseudo-spherical direct-beam correction. |
| 17 | a–b | blocked by solver feature | CPP has delta-M scaling, but not DISORT 4.0.99 delta-M-plus. |

The older cDISORT suite has 48 subcases.  Its Problems 1–13 overlap the
Fortran configurations above.  cDISORT Problem 14b additionally calls the
separate `twostr()` solver; CPP-DISORT with `NQuad=2` is not an implementation
of that algorithm and must not be presented as an exact port.
