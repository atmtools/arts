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
  take and give numbers.  Delta-M-plus scaling and the coupled calculations
  are DISORT's only, and Fresnel surfaces are VDISORT's only.  Both solvers
  have subsurface calculations, each with its own methods; that of VDISORT is
  described in `Subsurface emission with VDISORT`_.
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
* **Phase matrices from Legendre moments.**  VDISORT's settings hold the
  scattering of each layer as the Legendre moments of the six elements of the
  scattering matrix.  The phase matrices between the streams, and from the
  sun, are those of the scattering matrix's expansion in generalized
  spherical functions to the same degree.  This expansion agrees with the
  Legendre series in the phase function, :math:`F_{11}`, and differs from the
  truncated series of the other elements only by their truncation: it is
  regular at forward and backward scattering, where those truncated series
  are not.  See :doc:`dev.disort`.

Subsurface emission with VDISORT
--------------------------------

``spectral_radSubsurfaceVDisort`` gives the radiance that leaves the surface
at a point of a ray path from the emission and scattering of the subsurface
below it.  It is meant for :attr:`~pyarts3.workspace.Workspace.spectral_rad_surface_agenda`:

.. code-block:: python

   ws.subsurf_vdisort_settings_agendaSetup(scattering_setting="ScatteringSpecies")

   @pyarts.arts_agenda(ws=ws, fix=True)
   def spectral_rad_surface_agenda(ws):
       ws.spectral_radSubsurfaceVDisort(depth_profile=[0.0, -0.01, -0.1, -1.0])
       ws.spectral_radSurfaceDepolarization()

* **The subsurface.**  The temperature comes from
  :attr:`~pyarts3.workspace.Workspace.subsurf_field`, whose vertical
  coordinate is the altitude relative to the local surface (see
  :ref:`Sec Subsurface field`).  The absorption comes from the properties
  "refractive index", :math:`n'`, and "loss tangent", :math:`\tan\delta`, of
  :attr:`~pyarts3.workspace.Workspace.subsurf_spectral_field`, which give the
  complex refractive index :math:`n=\sqrt{n'^2(1+i\tan\delta)}`.  Volume
  scattering comes from
  :attr:`~pyarts3.workspace.Workspace.subsurf_scat_species`, evaluated at the
  points of the subsurface and its spectral field (see
  :doc:`user.scattering_species`).
* **The column.**  ``depth_profile`` gives the depths of the column below the
  ray path point, starting at 0 and descending, in metres.  The subsurface
  field must reach the deepest of them.  VDISORT's quadrature, Legendre and
  Fourier dimensions are those of the workspace.
* **The steps.**  :attr:`~pyarts3.workspace.Workspace.subsurf_vdisort_settings_agenda`
  sets up VDISORT for the column.  ``subsurf_vdisort_settings_agendaSetup``
  creates it from named options: scattering by the scattering species or
  none, the layers' thermal emission or none, blackbody radiation from below
  the deepest level or none, Fresnel reflection back into the subsurface at
  the interface or none, and an isotropic sky that enters through the
  interface, by default the cosmic background, or none.  The absorption is
  always included.
* **The interface.**  The surface at the ray path point is tilted by the
  slopes of the elevation of :attr:`~pyarts3.workspace.Workspace.surf_field`.
  The radiance within the subsurface that reaches the interface along the
  refracted direction is transmitted with one minus the Fresnel reflectance
  of each polarization.  To it is added the radiance of
  :attr:`~pyarts3.workspace.Workspace.spectral_rad_observer_agenda` along the
  direction of specular reflection, reflected by the Fresnel reflectance: the
  sky, with the sun where it is in that direction.
* **Depolarization.**  ``spectral_radSurfaceDepolarization`` scales the
  linear polarization of the outgoing radiance by :math:`1-D`, with
  :math:`D` the "depolarization factor" of
  :attr:`~pyarts3.workspace.Workspace.surf_spectral_field`.
* **Derivatives.**  The derivatives are those of the reflected sky by the
  targets of the observer agenda, and that of the depolarization factor.
  Derivatives by surface, subsurface and spectral subsurface targets are not
  computed by ``spectral_radSubsurfaceVDisort``, which refuses them.

Implementation and validation details are in :doc:`dev.disort`.
