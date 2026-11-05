.. _Sec Surface scattering user:

Surface scattering model
########################

ARTS represents the surface lower boundary of radiative transfer with a fully
polarimetric BRDF framework.  The physics — the boundary condition, the model
closed forms, Kirchhoff consistency, and the units convention — is developed
in :doc:`concept.surface_scattering`; the implementation is in
:doc:`dev.surface_scattering`.  This page describes the interface: what to
build, where to put it, and how to wire the agendas.

The three ingredients
*********************

Using the framework requires three things:

1. **Models** — instances of the surface scattering model classes, collected
   in the workspace variable :attr:`~pyarts3.workspace.Workspace.surface_models`
   (a :class:`~pyarts3.arts.MapOfSurfaceScatteringModel`).
2. **Masks** — per-model surface weights stored in
   :attr:`~pyarts3.workspace.Workspace.surf_field` under
   :class:`~pyarts3.arts.SurfacePropertyTag` keys *equal to the model names*.
3. **Agendas** — the ``SurfaceScatteringModel`` options of
   ``spectral_rad_surface_agenda``, which chain the scattering methods and
   call ``spectral_rad_incoming_agenda`` and
   ``spectral_rad_closed_surface_agenda`` for the incoming radiance and the
   surface emission.

.. graphviz::

  digraph g {
    bgcolor="#00000000";
    rankdir = "TD";
    ratio = auto;
    node [ color = "#0271BB" fontcolor = "white" style = "filled,rounded" shape = "rectangle" ];
    "models" [ label = "surface_models" style = "filled" ];
    "named lookup" [ label = "surface_models[name]" shape = "ellipse" ];
    "weighting" [ label = "surface_models.weighting_option" shape = "ellipse" ];
    "masks" [ label = "surf_field[SurfacePropertyTag(name)]" shape = "ellipse" ];
    "model classes" [ label = "Lambertian / Fresnel / FlatScalar" ];
    "raw weight" [ label = "mask value at the hit point" ];
    "agenda" [ label = "spectral_rad_surface_agenda" style = "filled" ];
    "models" -> "named lookup" [ arrowhead = "none" ];
    "models" -> "weighting" [ arrowhead = "none" ];
    "named lookup" -> "model classes";
    "masks" -> "raw weight";
    "raw weight" -> "weighting";
    "weighting" -> "agenda";
  }

The scattering models
*********************

There are three model families.  Each model class holds its data in a single
:class:`~pyarts3.arts.SurfaceSpectralData` member, which may vary with
frequency, latitude, and longitude (see `Spectral and spatial data`_):

.. list-table::
   :header-rows: 1

   * - Class
     - Physics
     - Data
   * - :class:`~pyarts3.arts.LambertianSurfaceScatterer`
     - Diffuse, isotropic; reflectivity :math:`r(\nu)`
     - ``reflectivity``
   * - :class:`~pyarts3.arts.FresnelSurfaceScatterer`
     - Specular; real refractive index :math:`n_2(\nu)` under medium :math:`n_1`
     - ``refractive_index``, and ``n1`` (a :class:`float`, default 1.0)
   * - :class:`~pyarts3.arts.FlatScalarSurfaceScatterer`
     - Specular; angle-independent reflectance :math:`R(\nu)`
     - ``reflectivity``

At run time the data is evaluated at the surface hit point for every frequency
of the simulation ``freq_grid``.  Reflectivity-like values are expected in
[0, 1] and are clamped; refractive indices and ``n1`` must stay positive.

Spectral and spatial data
*************************

A :class:`~pyarts3.arts.SurfaceSpectralData` is data over frequency [Hz],
latitude [deg], and longitude [deg], in that order.  It is implicitly
constructible from three kinds of data:

:class:`float`
  Constant at all frequencies and everywhere on the surface.

:class:`~pyarts3.arts.SpectralGeodeticField3`
  Gridded data with the grids ``(frequency, latitude, longitude)`` —
  *frequency first* — and data of shape ``(nfreq, nlat, nlon)``.  The data is
  interpolated linearly in all three dimensions.  The latitude grid must lie in
  [-90, 90] and the longitude grid in [-180, 180); 180 itself is not allowed.
  Single-point latitude and longitude grids make the data constant in space,
  which is how a plain spectrum is stored.

A function ``f(frq, lat, lon) -> float``
  A :class:`~pyarts3.arts.NumericTernaryOperator`, evaluated at every
  frequency and hit point.

A :class:`~pyarts3.arts.SurfaceSpectralData` can be called as
``data(frq, lat, lon)`` to inspect the value a model will see.

Outside the grids of a :class:`~pyarts3.arts.SpectralGeodeticField3`, the
extrapolation is set per dimension and side by the attributes ``frq_low``,
``frq_upp``, ``lat_low``, ``lat_upp``, ``lon_low``, and ``lon_upp``
(:class:`~pyarts3.arts.InterpolationExtrapolation`): ``None`` (outside the grid
is a user error), ``Zero``, ``Nearest`` (the edge value), or ``Linear``
(unlimited linear extrapolation).  ``set_extrapolation(mode)`` sets all six at
once.  The defaults follow from the data when it is constructed:

- a dimension with a single grid point is ``Nearest``, so the data is constant
  along it;
- a dimension with several grid points is ``None``, so a query outside the
  stored range raises a ``Limit breached`` error.

Constant and functional data are defined everywhere and never extrapolate.

.. tip::

   The stored spectral resolution is decoupled from the simulation grid: a
   model built on a coarse 3-point spectrum works at any simulation
   frequency, subject to the extrapolation rule.

Building models and registering them
************************************

.. code-block:: python

  import numpy as np
  import pyarts3 as pa
  arts = pa.arts

  # A diffuse Lambertian surface with 30% reflectivity
  lambert = arts.LambertianSurfaceScatterer(0.3)

  # A specular water-like surface
  water = arts.FresnelSurfaceScatterer(2.7, 1.0)

  ws.surface_models = arts.MapOfSurfaceScatteringModel()
  ws.surface_models.add("land", lambert)
  ws.surface_models.add("water", water)

The map behaves like a :class:`dict`: ``ws.surface_models["land"]``,
``ws.surface_models["land"] = …``, ``"land" in ws.surface_models``,
``len(…)``, and iteration over ``(name, model)`` pairs all work.

Surface masks and weighting
***************************

Each model's raw weight at a surface hit point is the value of
``surf_field`` under the :class:`~pyarts3.arts.SurfacePropertyTag` whose string
equals the model's name in the map.  **The tag name and the map key must
match**; a model with no matching tag has weight 0 everywhere.  Any
:class:`~pyarts3.arts.SurfaceData` representation works — a constant, a
:class:`~pyarts3.arts.GeodeticField2`, or a Python callable of ``(lat, lon)``
(see :doc:`user.surf_field`):

.. code-block:: python

  ws.surf_field[arts.SurfacePropertyTag("land")]  = 1.0
  ws.surf_field[arts.SurfacePropertyTag("water")] = lambda lat, lon: abs(lat) / 90.0

Masks must be non-negative; negative values abort with a user error.

The raw masks become model weights through ``surface_models.weighting_option``
(:class:`~pyarts3.arts.MapOfSurfaceScatteringModel.Weighting`):

``Maximum`` (default)
  Winner takes all: the model(s) with the largest mask value get weight 1
  (split equally on ties), all others 0.  *All-zero masks give every model
  weight 1/N* — a full blend, not zero.

``Average``
  Mask-weighted average: the raw masks normalized to sum 1.  *All-zero masks
  give all-zero weights* — no scattering and no surface emission at that
  point.

.. warning::

   The all-zero mask behaves *differently* between the two options.  Under
   ``Average`` a zero mask switches the surface off entirely; under
   ``Maximum`` it blends all models equally.  If your surface vanishes or
   behaves strangely at points where the masks are zero, this is why.

Agenda wiring
*************

The full setup for a sun-lit surface — every line is load-bearing:

.. code-block:: python

  ws.spectral_rad_space_agendaSet(option="SunOrCosmicBackground")  # without this the sun contributes exactly nothing
  ws.spectral_rad_incoming_agendaSet(option="Emission")
  ws.spectral_rad_closed_surface_agendaSet(option="Blackbody")
  ws.ray_path_observer_agendaSetGeometric()  # drives the refraction-aware sun-path search
  ws.spectral_rad_surface_agendaSet(option="SurfaceScatteringModel")
  ws.spectral_rad_observer_agendaSet(option="Emission")

The ``SurfaceScatteringModel`` options of ``spectral_rad_surface_agenda``
chain the scattering methods for you:

.. list-table::
   :header-rows: 1

   * - Option
     - Chain
     - Use for
   * - ``SurfaceScatteringModel``
     - Init + Diffuse + Specular + DiffuseDirect + SpecularDirect, with the
       sun and emission de-duplicated
     - Sky **and** sun beams, each counted exactly once
   * - ``SurfaceScatteringModelDiffuseOnly``
     - Init + Diffuse + Specular (no sun beams)
     - Passive thermal-IR / microwave, no sun
   * - ``SurfaceScatteringModelDirectOnly``
     - Init + DiffuseDirect + SpecularDirect (sun beams only)
     - Cheap optical glint / solar reflection only

The non-``Direct`` methods additionally require the quadrature grids
:attr:`~pyarts3.workspace.Workspace.zen_grid`,
:attr:`~pyarts3.workspace.Workspace.az_grid` and their weights
:attr:`~pyarts3.workspace.Workspace.zen_grid_weights`,
:attr:`~pyarts3.workspace.Workspace.az_grid_weights`; the ``Direct`` methods
require :attr:`~pyarts3.workspace.Workspace.suns` and
``ray_path_observer_agenda``.  The legacy options ``Blackbody``,
``Transmission``, and ``SurfaceReflectance`` are unrelated to this framework —
see :doc:`user.surface_emissivity` for the older reflectance path.

The quadrature grids
********************

The weights are applied verbatim, pairwise per direction — they must carry the
projected solid angle (see :ref:`Sec Surface scattering concept`).  Two
working recipes:

.. code-block:: python

  # Plain solid-angle hemisphere: weights sin(za) dza daa, sums to ~2*pi sr
  za = np.arange(1, 180, 2.0)
  aa = np.arange(0, 360, 2.0)
  ws.zen_grid = arts.ZenGrid(za.tolist())
  ws.az_grid  = arts.AziGrid(aa.tolist())
  ws.zen_grid_weights = arts.ZenGrid((np.sin(np.radians(za)) * np.radians(2.0)).tolist())
  ws.az_grid_weights  = arts.AziGrid(np.full(len(aa), np.radians(2.0)).tolist())

Or a cosine-hemisphere scheme whose weights already include the
:math:`\cos\theta` factor and sum to 1, which makes the scattered term a
direct weighted average of the incoming radiance.

.. note::

   On tilted surfaces extend ``zen_grid`` beyond 90°: the horizon gate uses
   the *actual* surface normal, and directions made visible by the tilt have
   zenith angles above 90° in surface coordinates.  A grid capped at 90°
   silently drops them.

   A coarse quadrature cannot resolve the ~0.5° solar disc.  Sunlight reaches
   the surface through the quadrature only by accident; use the ``Direct``
   methods (or the full option) for the sun.

Hand-written chains
*******************

Inside a ``@pa.workspace.arts_agenda`` bound to ``spectral_rad_surface_agenda``
you may chain the five methods —
:meth:`~pyarts3.workspace.Workspace.spectral_radSurfaceScatteringInit`,
``…Diffuse``, ``…Specular``, ``…DiffuseDirect``, ``…SpecularDirect`` — freely,
but then you own the de-duplication yourself:

- Run ``Init`` exactly once, at the top of the chain.  It sizes and zeroes
  ``spectral_rad`` / ``spectral_rad_jac``; every other method only adds.
  Omitting it raises a descriptive error.
- Set ``exclude_suns = 1`` on ``Diffuse``/``Specular`` whenever a ``Direct``
  method is in the same chain; keep it 0 otherwise, or the sun vanishes from
  the hemispheric route.
- Set ``include_emission = 0`` on the ``Direct`` methods whenever their
  non-``Direct`` counterpart is in the same chain; otherwise the emission is
  added twice.

With these set correctly, the identity
``full = DiffuseOnly + DirectOnly − one emission term`` holds to machine
precision whenever no traced direction contains a sun — a useful regression
check.

.. seealso::

   ``examples/7-surface-scattering-properties/`` contains runnable checkerboard
   scenes (pure surface, with atmospheric LUTs, IR, and camera variants).

Gotchas
*******

.. admonition:: Debugging traps seen in practice

   - **Emission double count**: a hand-written chain pairing a non-``Direct``
     method with its ``Direct`` counterpart without ``include_emission = 0``
     doubles the surface emission.
   - **Sun double count**: the same chain without ``exclude_suns = 1``
     doubles the sun; the reverse mistake — ``exclude_suns = 1`` under
     ``DiffuseOnly`` — deletes the sun entirely.
   - **Zero masks** behave differently under ``Maximum`` (full blend) and
     ``Average`` (surface off); negative masks abort.
   - **Weight convention**: pure ``sin(za) dza daa`` weights reflect a
     uniform radiance field with effective reflectance :math:`2\pi r`, not
     :math:`r`.  The weights must carry the projected solid angle.
   - **``zen_grid`` capped at 90°** silently loses the tilted-surface
     horizon.
   - **No sun without ``SunOrCosmicBackground``**: the ``Direct`` methods
     trace the sun through ``spectral_rad_incoming_agenda``; if the space
     agenda of that chain does not know the sun, the beam radiance is zero.
   - **No Jacobian with respect to model parameters or masks** — see
     :ref:`Sec Surface scattering concept`.
   - **``max_stepsize`` dividing ``toa`` exactly** (default 1000 m) can land a
     geometric path point a few 1e-10 m above the TOA, which the atmospheric
     field rejects.  Nudge either value if you hit this.
   - **XML round-trip gaps**: ``weighting_option`` is not serialised; re-set
     it after loading from XML.  Model data given as a Python function cannot
     be written to XML.
