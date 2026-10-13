DISORT input and output conventions
===================================

See :doc:`concept.disort` for the scalar and polarized transfer equations.

Gridded fluxes and DFDT are evaluated at each layer's lower boundary; the
top boundary is omitted.  See :doc:`user.heating` for conversion to K/s,
frequency integration, and pressure-based heating diagnostics.

* **Transparent layers.**  Exact zero-thickness layers are not accepted.
  ``tau_arr`` is a strictly increasing cumulative optical-depth grid beginning
  above zero, so each represented layer has positive thickness.  A transparent
  layer should be removed before constructing the solver.  Reference
  surface-only tests use a negligible but nonzero optical thickness and albedo;
  this approximates, but does not claim, exact transparent-layer support.
* **No separate two-stream solver.**  Setting :math:`N_q=2` selects a
  two-stream discrete-ordinate instance rather than a distinct approximate
  two-stream algorithm.
* **Phase conventions.**  VDISORT assumes all Mueller matrices have consistent
  propagation directions and Stokes reference planes.  The supplied combined
  matrices already include the cosine/sine transformation.  Applying that
  transformation twice, adding an extra factor of two, or mixing
  :math:`[I,Q,U,V]` with :math:`[I_v,I_h,U,V]` changes the physics.

Workspace methods
-----------------

* **The same methods for both solvers.**  DISORT's workspace methods and
  agendas have VDISORT counterparts named with ``vdisort`` in place of
  ``disort``.  Those of VDISORT take and give Stokes vectors where DISORT's
  take and give numbers.  Delta-M-plus scaling and the subsurface and
  coupled calculations are DISORT's only, and Fresnel surfaces are VDISORT's
  only.
* **Shared variables and options.**  Both solvers use the same quadrature,
  Legendre and Fourier dimensions, write the same quadrature and spectral
  flux field, and take the same named setup options.  VDISORT refuses the
  delta-M-plus scattering option.
* **Directions are propagation directions.**  The radiance fields are
  indexed by the direction in which the radiation propagates, as is the line
  of sight of a ray path point: a zenith angle below 90 degrees is upwelling
  radiation.  VDISORT's azimuth grid is the propagation azimuth too.
* **The sun.**  The settings' solar zenith and azimuth angles are those of
  the line of sight toward the sun.
* **Surface emission.**  Surface emission by temperature is the Planck
  function of the surface temperature, whatever the surface reflects.

Implementation and validation details are in :doc:`dev.disort`.
