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
  Legendre and Fourier dimensions, write the same quadrature and profiles of
  the flux components, and take the same named setup options.  VDISORT
  refuses the delta-M-plus scattering option.
* **Directions are propagation directions.**  The radiance fields are
  indexed by the direction in which the radiation propagates, as is the line
  of sight of a ray path point: a zenith angle below 90 degrees is upwelling
  radiation.  VDISORT's azimuth grid is the propagation azimuth too.
* **The sun.**  The settings' solar zenith and azimuth angles are those of
  the line of sight toward the sun.
* **Surface reflection.**  The surface methods (Lambertian, Hapke, Cox-Munk,
  RPV, Ross-Li, and Fresnel for VDISORT, each also spectral) store the
  surface's reflection between the streams, and its reflection of the sun's
  beam at the sun's zenith angle.  Set the sun before the surface: a solver
  with a beam refuses a surface computed for another sun.  A Fresnel surface
  cannot reflect a beam, so VDISORT refuses it with a sun.  For radiances at
  user directions, the low-level solvers also need the surface's reflection
  into those directions, for the same sun.
* **Surface emission.**  The upward boundary is the radiance beneath the
  surface: surface emission by temperature sets it to the Planck radiance
  of the surface temperature.  The solvers transmit it through the surface
  with one minus its reflection, which is Kirchhoff's emission: unpolarized
  :math:`(1-A)B` over a Lambertian surface and the polarized Fresnel
  emission over a flat one.  Set the surface reflection and the emission in
  either order.  Without a surface the boundary is the upward radiance
  itself.  See :doc:`concept.disort`.

Implementation and validation details are in :doc:`dev.disort`.
