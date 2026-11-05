.. _Sec Surface scattering dev:

Surface scattering framework development
########################################

This page describes the implementation of the surface scattering framework:
where the code lives, the contract a model must satisfy, how to add a new
model, and how the methods and agendas are assembled.  The physics is in
:doc:`concept.surface_scattering`, the interface in
:doc:`user.surface_scattering`.

Where the code lives
********************

.. list-table::
   :header-rows: 1

   * - Concern
     - File(s)
   * - Model bodies
     - ``src/core/surface_scattering/{lambertian,fresnel,flat_scalar}.{h,cc}``
   * - Model data over (frequency, latitude, longitude),
       ``SurfaceSpectralData``
     - ``src/core/surface_scattering/surface_spectral_field.{h,cc}``
   * - Field interpolation and per-dimension extrapolation
     - ``src/core/functional/functional_field_interp.h``
   * - Variant type, named map, weighting, blending
     - ``src/core/surface_scattering/surface_scattering_model.{h,cc}``
   * - Output tensor types
     - ``src/core/surface_scattering/surface_scattering_properties.h``
   * - Workspace method implementations
     - ``src/m_surface.cc``
   * - Method registry (in/out/gin declarations)
     - ``src/workspace_methods.cpp``
   * - Agenda contracts and option enums
     - ``src/workspace_agendas.cpp``
   * - Predefined option → method chains and gin presets
     - ``src/workspace_agenda_creator.cpp``
   * - Python bindings (the three models, ``SurfaceSpectralData``, and the
       map)
     - ``src/python_interface/py_surface_scattering.cpp``
   * - Tests
     - ``tests/core/surf/`` (collected by ``collect_test_subdir()``)

The source files are enumerated in ``src/core/surface_scattering/CMakeLists.txt``
(not globbed) — a new model source must be added there.

The properties contract
***********************

Every model answers one question — *what are your scattering properties at
this point, for these frequency and angle grids?* — via a const member with a
fixed signature:

.. code-block:: cpp

   SurfaceScatteringModelProperties
   get_surface_scattering_model_properties(
       const SurfacePoint& surf_point,
       Numeric lat, Numeric lon,
       const Vector& f_grid,
       const Vector& za_inc_grid, const Vector& aa_inc_grid,
       const Vector& za_scat_grid, const Vector& aa_scat_grid) const;

The return type, ``surface_scattering::SurfaceScatteringModelProperties``,
holds four ``std::optional`` components:

.. list-table::
   :header-rows: 1

   * - Component
     - Type and dimensions
   * - ``brdf_matrix_diffuse`` / ``brdf_matrix_specular``
     - ``BRDFMatrix`` = ``MuelmatTensor5``,
       ``[nf, nza_inc, naa_inc, nza_scat, naa_scat, 4, 4]``
   * - ``emissivity_vector_diffuse`` / ``emissivity_vector_specular``
     - ``EmissivityVector`` = ``MuelmatTensor3``,
       ``[nf, nza_scat, naa_scat, 4, 4]``

A model fills only the pair of its own family and leaves the other pair
``std::nullopt``; consumers read a missing component as a zero contribution.
The struct supports ``operator+=`` (accumulating across models) and
``operator*=`` (applying the blend weight).

The map dispatches structurally, not by name
********************************************

``MapOfSurfaceScatteringModel`` visits the variant with a generic lambda
wrapped in a ``static_assert``:

.. code-block:: cpp

   static_assert(requires {
                   model.get_surface_scattering_model_properties(
                       surf_point, lat, lon, f_grid,
                       za_inc_grid, aa_inc_grid, za_scat_grid, aa_scat_grid);
                 },
                 "Every surface scattering model variant member must implement "
                 "get_surface_scattering_model_properties");

Any variant alternative that satisfies the contract compiles; any that does
not fails at compile time with a named error.  Consequently **adding a model
requires no edit** to ``m_surface.cc``, to the five workspace methods, or to
any registry file.

Adding a new scattering model
*****************************

1. Write the struct under ``src/core/surface_scattering/` with the member
   above, returning the tensor pair of your family at the requested shapes.
   Hold spectrally or spatially varying inputs as ``SurfaceSpectralData``
   members and evaluate them with their ``at(f_grid, lat, lon)``.  Add the
   source to ``CMakeLists.txt``.
2. Append the type to the ``SurfaceScatteringModel`` ``std::variant`` in
   ``surface_scattering_model.h`` and add a global ``using`` alias (mirroring
   the existing three).
3. Add a ``std::formatter`` specialisation — without it the model is invisible
   to Python ``print`` — and an ``xml_io_stream`` specialisation with a unique
   ``type_name``.  The variant XML machinery requires default-constructibility
   and enforces it at compile time.
4. Register a ``py::class_`` in ``py_surface_scattering.cpp`` and call
   ``generic_interface(…)`` on it — the nanobind variant caster converts only
   registered alternatives, so an unregistered model cannot be passed to
   ``surface_models.add``.
5. Rebuild.  The model is immediately usable from Python; add a per-model unit
   test under ``tests/core/surf/``.

Methods, agendas, and the de-duplication gin
********************************************

The five workspace methods are declared in ``src/workspace_methods.cpp`` and
implemented in ``src/m_surface.cc``:

.. list-table::
   :header-rows: 1

   * - Method
     - Samples
     - Gin (default)
   * - ``spectral_radSurfaceScatteringInit``
     - — (sizes and zeroes the outputs)
     - none
   * - ``spectral_radSurfaceScatteringDiffuse``
     - ``zen_grid`` × ``az_grid`` quadrature
     - ``exclude_suns`` (0)
   * - ``spectral_radSurfaceScatteringSpecular``
     - single mirror direction (no quadrature)
     - ``exclude_suns`` (0)
   * - ``spectral_radSurfaceScatteringDiffuseDirect``
     - one delta beam per sun
     - ``angle_cut`` (0.0), ``refinement`` (1), ``include_emission`` (1)
   * - ``spectral_radSurfaceScatteringSpecularDirect``
     - sun disc test on the mirror direction
     - ``angle_cut`` (0.0), ``refinement`` (1), ``include_emission`` (1)

All scattering methods read ``freq_grid``, ``atm_field``, ``surf_field``,
``subsurf_field``, ``surface_models``, ``jac_targets``, ``ray_point``,
``suns``, and both sub-agendas; the ``Direct`` methods additionally require
``ray_path_observer_agenda`` for the refraction-aware sun-path search (the
same iterative solver as ``sun_pathFromObserverAgenda`` with ``just_hit = 1``;
``angle_cut`` and ``refinement`` are its control parameters).

The predefined agenda options are assembled in
``get_spectral_rad_surface_agenda`` in ``src/workspace_agenda_creator.cpp``:

- ``SurfaceScatteringModel``: Init, then Diffuse and Specular with
  ``SetWsv("exclude_suns", Index{1})``, then DiffuseDirect and SpecularDirect
  with ``SetWsv("include_emission", Index{0})``.
- ``SurfaceScatteringModelDiffuseOnly``: Init + Diffuse + Specular at gin
  defaults.
- ``SurfaceScatteringModelDirectOnly``: Init + DiffuseDirect + SpecularDirect
  at gin defaults.

``AgendaCreator::add(name, SetWsv(...))`` encodes a gin preset as a
``@``-prefixed named-input entry in the finalized agenda (``named_input_prefix
= '@'``) — filter ``@``- and ``_``-prefixed method names when comparing chains
name-by-name.

Tests
*****

Run with ``ctest --test-dir build -R 'tests.core.surf'``.  The layout under
``tests/core/surf/``:

- **Model unit tests** call ``get_surface_scattering_model_properties``
  directly on a model or the map and assert tensor shapes and closed-form
  values: ``lambertian_scatterer_field.py``, ``fresnel_surface_scatterer.py``,
  ``flat_scalar_surface_scatterer.py``, ``interp_extrapolation.py``,
  ``map_weighting_option.py``.  This is the cheapest way to test a new model —
  no workspace, no agendas.
- **Per-method tests** with closed-form radiance checks:
  ``spectral_rad_surface_scattering_{diffuse,specular,diffuse_direct,specular_direct}.py``.
- **Chaining and option tests**: ``spectral_rad_surface_scattering_chain.py``,
  ``spectral_rad_surface_scattering_predefined_agendas.py``,
  ``spectral_rad_surface_scattering_sun_double_count.py``.
- **Sub-surface coupling**: ``spectral_radiance_subsurface.py`` and the
  ``spectral_radiance_jacobian_subsurface*.py`` family.

Known serialization gaps
************************

- ``MapOfSurfaceScatteringModel`` XML I/O writes only the models;
  ``weighting_option`` is not serialised.

This is a candidate fix rather than intended behavior — extend the
corresponding ``xml_io_stream`` specialisation if you depend on it.
