Radiative heating rates
=======================

``pyarts3.recipe.heating_rates`` converts radiance and flux profiles to gas
temperature tendencies in K/s. Positive values mean warming; multiply by
86400 for K/day. Supply mass specific heat capacity in J/(kg K), density in
kg/m³ where required, pressure in Pa, and gravity in m/s². The recipe accepts
ARTS arrays and NumPy arrays; it does not infer units from their values.

From radiance to flux
---------------------

Heating requires net **projected** irradiance, not angle-integrated radiance
or actinic flux. ``flux_from_radiance`` integrates Stokes I with the zenith
cosine and supplied quadrature weights. Input radiance must be in
W/(m² sr Hz), or W/(m² sr) after frequency integration. Brightness temperature
must first be converted to physical radiance. The helper supports azimuth
independent radiation or an azimuth average, and uses ARTS viewing zenith
angles: a ray looking upward measures downward travelling radiation.

For an array with axes (frequency, altitude, zenith)::

    from pyarts3.recipe import heating_rates as heating

    up, down = heating.flux_from_radiance(
        spectral_radiance_I, zenith, zenith_weights)
    net = heating.integrate_spectral(up - down, frequency)
    rate = heating.from_flux(net, pressure, heat_capacity, gravity)

``zenith_weights`` integrate over cosine, for example the double-Gauss
weights. Frequency integration uses trapezoids on an increasing Hz grid or
explicit quadrature weights in Hz. Frequency and angle integration can be
exchanged. A direct solar beam is a separate irradiance contribution and must
be added to the downward flux if absent from the diffuse radiance field.
For an ARTS gridded Stokes field, select Stokes I and arrange the axes as above
(or pass the appropriate ``axis``); no global ARTS2 spatial grids are needed.

From flux to heating
---------------------

Under hydrostatic balance, with upward net flux :math:`F=F_\uparrow-F_\downarrow`,

.. math::

   \frac{dT}{dt} = \frac{g}{c_p}\frac{dF}{dp}.

``from_flux`` evaluates a three-point derivative on nonuniform pressure
levels, with second-order one-sided boundaries. Pressure can increase or
decrease. Gravity and heat capacity can vary with location. For a different
chosen pressure-differencing scheme, ``from_flux_divergence`` converts an
already computed :math:`dF/dp`. This plane-parallel vertical diagnostic does
not include horizontal flux divergence or spherical-area divergence.

From DISORT optical-depth derivatives
-------------------------------------

DISORT and VDISORT return DFDT with respect to **optical depth**, not time:

.. math::

   D_\nu = \frac{dF_\nu}{d\tau_\nu},\qquad
   \frac{dT}{dt} = \frac{1}{\rho c_p}\int k_{\mathrm{ext},\nu}D_\nu\,d\nu.

Total extinction is required: the absorption fraction is already present in
DFDT. Frequency-dependent extinction must be applied before integration::

    rate = heating.from_disort(
        ws.spectral_flux_components_profile,
        extinction, density, heat_capacity,
        weights=frequency_quadrature_weights)

Omit ``weights`` for trapezoidal integration. Extinction broadcasts to
(frequency, layer); density and heat capacity normally have shape (layer,).
Workspace ``FluxComponentsProfile`` stores values at ``alt_grid[1:]``, the
lower boundary of each layer, and uses optical properties of the layer
immediately above that boundary. It omits the top boundary. It does not store
layer averages or layer-center values. Thermodynamic inputs must match the
output locations. This local DFDT diagnostic differs from finite-differencing
flux on pressure levels and need not produce identical numbers.

For low-level ``cppdisort`` or ``cppvdisort`` output, select row 3 of
``solver.flux(tau)`` and call ``from_optical_depth_derivative``. This preserves
spectral dimensions; integrate the resulting K/(s Hz) values afterwards.
Extinction must correspond to the optical-depth coordinate of DFDT. Scalar
DISORT's internal delta-M treatment returns physical DFDT, so use physical
extinction. If VDISORT receives externally rescaled transport inputs, use the
extinction corresponding to those inputs rather than mixing coordinates.


See :doc:`concept.heating` for the energy-balance equations and
:doc:`dev.heating` for reference-test adaptations.
