Workspace methods
#################

Workspace methods are the core interaction a user has with ARTS calculations.
They are called from the :class:`~pyarts3.workspace.Workspace` object in python but are
implemented in C++ under-the-hood.

The workspace methods definition are mainly located in the ``workspace_methods.cpp``
and ``workspace_meta_methods.cpp`` files.  Some workspace methods are 
automatically generated elsewhere.

Workspace methods are exposed to other C++ files via the ``<workspace.h>`` header file.
Only the CMake target ``artsworkspace`` and targets that depend on it may include
``<workspace.h>``.

How to add a workspace method?
==============================

The easiest way to add a workspace method is to copy an existing method.
Go to the ``workspace_methods.cpp`` or ``workspace_meta_methods.cpp`` file, depending on
if it is a basic or meta-method you wish to add.
Copy an existing method and modify it to suit your needs.

Meta-methods are automatically generated and do not require any additional code.
Basic methods require new code to be compiled as part of the ``artsworkspace`` target.
Please add this code to an existing ``m_``-file or create a new one.
Please ensure that the new file is listed in the CMake target ``artsworkspace``.

Defining a workspace method
===========================

You may define either a basic method or a meta-method that depends on other methods.
The easiest way to do so is to copy an existing method and modify it to suit your needs.
Please design your method input and output to allow meta-methods to be defined,
as they simplify composability significantly.

Basic methods
-------------

Basic methods completely define their input and output and do not explicitly
depend on other methods.  The method must be manually implemented.

Basic methods are defined entirely in ``workspace_methods.cpp``.
They are part of a map object called ``wsm_data``.  The name of the
method is the key and the value is a struct with the following fields:

- ``desc`` - a description of the method as a string.
- ``author`` - the author(s) of the method as a list of strings.
- ``out`` - the :doc:`dev.workspace.variables` output of the method as a list of strings.
- ``gout`` - the generic output of the method as a list of strings.  These must not be :doc:`dev.workspace.variables`.
- ``gout_type`` - the type of the generic output as a list of strings.  These must be :doc:`dev.workspace.groups`.
- ``gout_desc`` - the description of the generic output as a list of strings.
- ``in`` - the :doc:`dev.workspace.variables` input of the method as a list of strings.
- ``gin`` - the generic input of the method as a list of strings.  These must not be :doc:`dev.workspace.variables`.
- ``gin_type`` - the type of the generic input as a list of strings.  These must be :doc:`dev.workspace.groups`.
- ``gin_value`` - the default value of the generic input as an optional initialized :doc:`dev.workspace.groups`.
- ``gin_desc`` - the description of the generic input as a list of strings.
- ``pass_workspace`` - a boolean indicating if a :class:`~pyarts3.workspace.Workspace` instance should be passed to the method.  If true, the first argument to the method is a ``const Workspace&``.

The expected signature of the method will depend on these fields.
A linker error will likely occur if the actual signature does not match
the expected signature.

The ``in`` and ``out`` may contain the same :doc:`dev.workspace.variables`.  If they do, the variable must be
initialized before the method is called because it is treated as if the method is intended to
simply modify the existing value.  Please indicate strongly in the documentation if you sometimes overwrite the input variable.

On the other hand, if the :doc:`dev.workspace.variables` is only in ``out`` and not in ``in``,
it is treated as if the workspace variable is created by the method.  Note that since the type system
does not account for this, it is important that you clear the current state of the :doc:`dev.workspace.variables`
in a method that is intended to create a new workspace variable.

The fields ``gin``, ``gin_type``, ``gin_value``, and ``gin_desc`` must be the same size.
The same is true for ``gout``, ``gout_type``, and ``gout_desc``.  These are user-generated
inputs and outputs, and are often used to pass information pertinent to the method itself
but not to the workspace as a whole.

Please check other workspace methods for examples by comparing their actual signature
to the expected signature to figure out how the fields should be filled in.  Also check
that the documentation is generated as intended by building the ``pyarts_docs`` target.

.. tip::

  All fields but ``desc`` and ``author`` are optional.  If a field is not needed, it
  is convenient to leave it out.

Meta-methods
------------

Meta-methods do not define all their input and output, but instead define a call
order into other methods.  From this call order, the inputs of the user-facing
workspace method is inferred.  This method should not be implemented manually.

These methods are defined in ``workspace_meta_methods.cpp``.  They are defined
as part of a list called ``wsm_meta``.
A single meta-method data contains:

- ``name`` - the name of the method as a string.
- ``desc`` - a description of the method as a string.
- ``author`` - the author(s) of the method as a list of strings.
- ``methods`` - the methods that the meta-method depends on as a list of strings.
- ``out`` - the output of the method as a list of strings.  These must be workspace variables.
- ``preset_gin`` - The preset ``gin`` values for the method as a list of workspace values.
- ``preset_gin_value`` - The preset ``gin_value`` values for the method as a list of workspace values.

.. tip::

  A meta-method may depend on another meta-method.  If it does, it is important that the
  meta-method it depends on is defined before it in the list.

Automatic methods
-----------------

All methods that execute a workspace agenda are automatically generated.
These will be named as ``agenda_nameExecute`` and may otherwise be
treated as normal workspace method.
You need to do nothing to define these methods.  But please refrain from defining
them manually as that may cause undefined naming conflicts.

The expected signature of the method :func:`~pyarts3.workspace.Workspace.spectral_propmat_agendaAuto` is also
generated automatically near the end of ``workspace_methods.cpp``.  It takes
its input and output from a list of other methods.  Feel free to add to this
list but make sure that any naming conflicts regarding ``gin`` are resolved
before doing so.  Adding a method to this list may also require changing the
actual signature (which is why the method is generated, so that a change in
the required actual signature is immediately made apparent).

The methods that begin with ``RetrievalAdd...`` are partly generated.
These methods all require a corresponding ``jac_targetsAdd...`` method
that fills in the ``jac_targets`` workspace variable.  To keep that
part of the signature consistent, the additional ``RetrievalAdd...`` information
is simply appended to the ``in``, ``out``, and ``gin``-lists of the
corresponding ``jac_targetsAdd...`` method using the local ``jac2ret`` lambda.

Generated files
===============

The workspace method interface generates a lot of files during the build process.
These generated files are located in the build directory and are named
as ``auto_wsm_N.cc``, where N is a number, as ``auto_wsm.cpp``, as ``auto_wsm.h``,
and as ``auto_wsmmeta.cpp`` for the C++ interfacing code.  The python-binding
code is also generated as ``py_auto_wsm_N.cpp``, where N is still a number.

Workspace method naming convention
==================================

Names carry meaning.  Please follow the naming convention below, and
please do not hesitate to fix any naming inconsistencies you find.

Method naming
-------------

Workspace method names should be descriptive and follow the naming convention
that the main workspace variable output of the method in ``snake_case``
is followed by a short but descriptive name of what the method does with the output
in ``PascalCase``.
A general rule of thumb is to use verbs for methods that modify the workspace
variable and nouns for methods that create a new workspace variable.

For example, :func:`~pyarts3.workspace.Workspace.spectral_propmatAddLines`
has a main output of :attr:`~pyarts3.workspace.Workspace.spectral_propmat` and
adds line absorption to it.  It needs to be preceded by a call to 
:func:`~pyarts3.workspace.Workspace.spectral_propmatInit` which sets up the
propagation matrix to an initial state.

Of course, every use-case is different, but please try to follow this convention.

File naming
-----------

The file that a workspace method is implemented in should be named ``m_<concept>.cc``.
The concept should be a short but descriptive name of what the methods therein do.
Multiple methods per file is allowed and encouraged, but keep them conceptually similar.
To ensure compatibility with various file systems, please avoid using spaces
and capital letters in the filename.

Lastly, please ensure that the file is listed in the CMake target ``artsworkspace``,
or it will not be compiled.

Workspace method documentation
==============================

Workspace documentation that contains ``*text*`` is automatically turned into links
to the relevant ARTS-related variable or method.  Please use this feature to link
between workspace methods and variables.

If a method require extra information beyond what you can fit in the ``desc`` field,
there's a ``workspace_method_extra_doc.cpp`` file that you can add to.  This file
has access to the full workspace as part of the ``artsworkspace`` target and the 
python documentation adds a separate subsection for the information in this file (documentation level ``-------``).

Examples of defined workspace methods
=====================================

The following examples are taken from the ARTS source code.  Please check the
source code for the full context of the examples.

Method creating a workspace variable
------------------------------------

The following is a basic
method that creates or set a workspace variable.

This is the extration of the text in the ``workspace_methods.cpp`` file:

.. code-block:: c++

    wsm_data["ray_pathGeometricUplooking"] = {
        .desc =
            R"--(Wraps *ray_pathGeometric* for straight uplooking paths from the surface altitude at the position
    )--",
        .author = {"Richard Larsson"},
        .out    = {"ray_path"},
        .in     = {"atm_field", "surf_field", "lat", "lon", "max_stepsize"},
    };

The signature of the method is:

.. code-block:: c++

  void ray_pathGeometricUplooking(ArrayOfPropagationPathPoint& ray_path,
                                  const AtmField& atm_field,
                                  const SurfaceField& surf_field,
                                  const Numeric& lat,
                                  const Numeric& lon,
                                  const Numeric& max_stepsize);

The signature of the method returns ``void``.  This is the same for all ARTS methods.

The first argument of the method is a reference to :attr:`~pyarts3.workspace.Workspace.ray_path`.
Since :attr:`~pyarts3.workspace.Workspace.ray_path` is in ``out`` but not in ``in``,
it is expected that the method overwrite any existing value of :attr:`~pyarts3.workspace.Workspace.ray_path`.

The arguments :attr:`~pyarts3.workspace.Workspace.atm_field`, :attr:`~pyarts3.workspace.Workspace.surf_field`,
:attr:`~pyarts3.workspace.Workspace.latitude`, and :attr:`~pyarts3.workspace.Workspace.longitude`
are defined in ``in`` and are passed to the method as immutable references to the respective
workspace variables.

All other fields are there to provide context and to generate the documentation.
See :meth:`~pyarts3.workspace.Workspace.ray_pathGeometricUplooking` for the full documentation.

.. note::
  The command ``*ray_pathGeometric*``
  is automatically turned to a link by the pyarts build system
  by mimicying of the sphinx linking system.

Method modifying a workspace variable
-------------------------------------

The following is a basic workspace method that modifies existing workspace variables.

This is the extraction of the text in the ``workspace_methods.cpp`` file:

.. code-block:: c++

    wsm_data["spectral_propmatAddLines"] = {
        .desc      = R"--(Add line-by-line absorption to the propagation matrix.

    See :doc:`concept.absorption.lbl` for details.
    )--",
        .author    = {"Richard Larsson"},
        .out       = {"spectral_propmat",
                      "spectral_nlte_srcvec",
                      "spectral_propmat_jac",
                      "spectral_nlte_srcvec_jac"},
        .in        = {"spectral_propmat",
                      "spectral_nlte_srcvec",
                      "spectral_propmat_jac",
                      "spectral_nlte_srcvec_jac",
                      "freq_grid",
                      "jac_targets",
                      "select_species",
                      "abs_bands",
                      "abs_ecs_data",
                      "atm_point",
                      "ray_point"},
        .gin       = {"no_negative_absorption"},
        .gin_type  = {"Index"},
        .gin_value = {Index{1}},
        .gin_desc =
            {"Turn off to allow individual absorbers to have negative absorption"},
    };

The signature of the method is:

.. code-block:: c++

  void spectral_propmatAddLines(PropmatVector& spectral_propmat,
                                StokvecVector& spectral_nlte_srcvec,
                                PropmatMatrix& spectral_propmat_jac,
                                StokvecMatrix& spectral_nlte_srcvec_jac,
                                const AscendingGrid& freq_grid,
                                const JacobianTargets& jac_targets,
                                const SpeciesEnum& select_species,
                                const AbsorptionBands& abs_bands,
                                const LinemixingEcsData& abs_ecs_data,
                                const AtmPoint& atm_point,
                                const PropagationPathPoint& ray_point,
                                const Index& no_negative_absorption);

The signature of the method returns ``void``.  This is the same for all ARTS methods.

The first four arguments of the method are references to
:attr:`~pyarts3.workspace.Workspace.spectral_propmat`.
:attr:`~pyarts3.workspace.Workspace.spectral_nlte_srcvec`,
:attr:`~pyarts3.workspace.Workspace.spectral_propmat_jac`, and
:attr:`~pyarts3.workspace.Workspace.spectral_nlte_srcvec_jac`
are both output (``out``) and input (``in``).  The method is expected to modify the existing values
of these workspace variables instead of creating new ones.

The arguments
:attr:`~pyarts3.workspace.Workspace.freq_grid`,
:attr:`~pyarts3.workspace.Workspace.jac_targets`,
:attr:`~pyarts3.workspace.Workspace.select_species`,
:attr:`~pyarts3.workspace.Workspace.absorption_bands`,
:attr:`~pyarts3.workspace.Workspace.abs_ecs_data`,
:attr:`~pyarts3.workspace.Workspace.atm_point`, and
:attr:`~pyarts3.workspace.Workspace.ray_point` are just defined in ``in`` and are passed to the method
as immutable references to the respective workspace variables.

Lastly, the argument ``no_negative_absorption`` is defined in ``gin`` and is passed
as an immutable reference as well.  The type of the argument is ``Index``
and the default value is ``1``.  The default value is passed to the method
if the user does not provide a value for ``no_negative_absorption``.
The ``no_negative_absorption`` argument is used to turn off the check for negative absorption,
which is useful for debugging purposes.

The other fields are there to provide context and to generate the documentation.
See :meth:`~pyarts3.workspace.Workspace.spectral_propmatAddLines` for the full documentation.

.. note::
  The link using the ``:doc:``-command is automatically generated by
  the sphinx documentation system.  Any default sphinx command, and
  a few extra commands beyond that, are automatically respected in
  full.  The only exception is the ``*text*``-command, which is
  automatically turned into links to ARTS-related variables and methods.

Method that uses a workspace agenda
-----------------------------------

The following is a basic workspace method that creates workspace variables.

This is the extraction of the text in the ``workspace_methods.cpp`` file:

.. code-block:: c++

  wsm_data["measurement_vecFromSensor"] = {
        .desc =
            R"--(Sets measurement vector by looping over all sensor elements

  The core calculations happens inside the *spectral_rad_observer_agenda*.

  User choices of *spectral_rad_unit* does not adversely affect this method
  unless the *measurement_vec* or *measurement_jac* are further modified
  before consumption by, e.g., *OEM*
  )--",
        .author         = {"Richard Larsson"},
        .out            = {"measurement_vec", "measurement_jac"},
        .in             = {"measurement_sensor",
                           "jac_targets",
                           "atm_field",
                           "surf_field",
                           "spectral_rad_unit",
                           "spectral_rad_observer_agenda"},
        .pass_workspace = true,
    };

The signature of the method is:

.. code-block:: c++

  void measurement_vecFromSensor(const Workspace& ws,
                                 Vector& measurement_vec,
                                 Matrix& measurement_jac,
                                 const ArrayOfSensorObsel& measurement_sensor,
                                 const JacobianTargets& jac_targets,
                                 const AtmField& atm_field,
                                 const SurfaceField& surf_field,
                                 const SpectralRadianceUnitType& spectral_rad_unit,
                                 const Agenda& spectral_rad_observer_agenda);

The signature of the method returns ``void``.  This is the same for all ARTS methods.

The first argument of the method is a reference to the workspace object itself.
This is passed as a ``const Workspace&`` reference to the method.  It is passed
to the method because ``pass_workspace`` is set to ``true`` in the method definition.
Note that the workspace object is passed as a ``const`` reference, so it cannot be modified.

The coming two arguments of the method are references to
:attr:`~pyarts3.workspace.Workspace.measurement_vec` and
:attr:`~pyarts3.workspace.Workspace.measurement_jac`.
Since :attr:`~pyarts3.workspace.Workspace.measurement_vec` and
:attr:`~pyarts3.workspace.Workspace.measurement_jac` are in ``out`` but not in ``in``,
it is expected that the method overwrite any existing values they might hold.

The arguments :attr:`~pyarts3.workspace.Workspace.measurement_sensor`,
:attr:`~pyarts3.workspace.Workspace.jac_targets`,
:attr:`~pyarts3.workspace.Workspace.atm_field`,
:attr:`~pyarts3.workspace.Workspace.surf_field`,
:attr:`~pyarts3.workspace.Workspace.spectral_rad_unit`, and
:attr:`~pyarts3.workspace.Workspace.spectral_rad_observer_agenda`
are defined in ``in`` and are passed to the method
as immutable references to the respective workspace variables.

The other fields are there to provide context and to generate the documentation.
See :meth:`~pyarts3.workspace.Workspace.measurement_vecFromSensor` for the full documentation.

Meta-method output with workspace variables
-------------------------------------------

The following is a meta-method that creates workspace variables.

This is the extraction of the text in the ``workspace_meta_methods.cpp`` file:

.. code-block:: c++

  wsm_meta.push_back(WorkspaceMethodInternalMetaRecord{
      .name             = "atm_fieldRead",
      .desc             = "Reads absorption file from a directory",
      .author           = {"Richard Larsson"},
      .methods          = {"atm_fieldInit",
                           "atm_fieldAppendBaseData",
                           "atm_fieldAppendAuto"},
      .out              = {"atm_field"},
      .preset_gin       = {"replace_existing"},
      .preset_gin_value = {Index{0}},
  });

The signature of the generated meta-method is
not important because it is generated automatically.

Calling the above method is effectively the same as calling
the listed methods one after the other and then deleting all method output
that is not in ``out``.
In other words, even if a sub-method has an output that is not in ``out``,
it will not be passed to the user.

The call order and documentation is available here
see :meth:`~pyarts3.workspace.Workspace.atm_fieldRead` 

Generic argument variants
-------------------------

A comma-separated generic type declaration lists independent alternatives
for that argument.  The generator produces one C++ signature using
``const Generic<const T, const U, ...>`` for inputs and
``Generic<T, U, ...>`` for outputs, passed by value.  ``Generic<Ts...>``
derives from ``std::variant<std::shared_ptr<Ts>...>``.  The outer ``const``
protects the input variant; the const-qualified alternatives protect its
pointees.  The ``UniformGenericConstness`` concept requires all alternatives to be
const or all to be non-const. Mixed alternatives are rejected.
``Any`` expands to ``AnyInput`` or ``AnyOutput``, covering all workspace
groups.  Workspace adapters retain the original shared pointers; they do
not copy the underlying values.  Inputs are const and outputs retain the
identity of the workspace object.  Assigning a different pointer to an
output variant does not replace the workspace variable: mutate its pointee.

For example, two ``AtmKey,SpeciesEnum`` inputs give two two-alternative
variants in one implementation, accepting all four combinations.  Visit
each input independently when converting to a common representation.
Use a multi-variant visit only when the algorithm actually depends on the
combination.  Cross-argument restrictions must be checked explicitly in
the implementation; type lists are no longer zipped overload tables.

The Python signature uses untyped objects for generic arguments.  One
shared runtime conversion implementation preserves explicitly typed
objects, then tries permitted implicit input conversions in declaration
order.  Generic outputs require an explicitly typed mutable ARTS object;
converting a Python scalar would otherwise mutate a discarded temporary.
No nanobind variant caster is instantiated per generic method signature.
The C++ adapter constructs the correctly sized variant from the resulting
workspace value.  Converted inputs retain ownership through the call.

For direct C++ calls, constructors accepting references borrow those
objects for the call's duration.  Such borrowed alternatives must not
escape the call.  Construct from workspace values or owning shared pointers
when retaining an argument beyond that lifetime.

Value variants also construct a ``Generic``, with deduction of their alternatives::

    Atm::KeyVal key = AtmKey::t;
    Generic mutable_key(key);
    Generic input_key(std::as_const(key));

Lvalue variants borrow their active value without copying it.  The source
must remain alive and keep the same active alternative while the borrowed
Generic is used.  Rvalue variants instead allocate shared ownership of the
moved active value; const rvalues copy it into owned const storage.
Const source variants produce const alternatives.

Combine alternative sets with ``Extend``::

    using First = Generic<AtmKey, SpeciesEnum>;
    using Second = Generic<SpeciesEnum, SpeciesIsotope>;
    using Combined = First::Extend<Second>;
    // Combined is Generic<AtmKey, SpeciesEnum, SpeciesIsotope>.
    AtmKey temperature = AtmKey::t;
    First first(temperature);
    Combined combined(first);

``Extend`` accepts multiple Generic types, preserves declaration order, and
removes exact duplicate types. All combined alternatives must have the
same constness; use ``Const`` before combining mutable and const types.
Combination extends the set of permitted types; the result still holds
one active value.  Converting an existing Generic copies or moves its
shared pointer, preserving ownership and object identity.  Conversions
can add pointee constness but cannot remove it.  All source alternatives
must be supported by the destination, including inactive alternatives.

Use ``Generic<T, U>::Const`` to obtain ``Generic<const T, const U>`` and
``value.as_const()`` to construct that type without copying its pointee.
An lvalue shares ownership; ``std::move(value).as_const()`` transfers its
active shared pointer.  Already const alternatives stay const.  This makes the pointees
read-only through the returned Generic; mutable aliases remain usable.
Merely declaring the outer Generic ``const`` does not make its pointees const.

``test_workspace_variant.cc`` checks type sizes, constness, ownership and
identity.  ``tests/core/agenda/supergeneric.py`` checks Python conversion
and mutable ``Any`` outputs, including execution through an agenda.

Keeping key alternatives consistent
-----------------------------------

Use a named helper such as ``AtmKeyValStr()`` in
``src/workspace_methods.cpp`` instead of repeating atmospheric key lists
in method metadata::

    .gin_type = {AtmKeyValStr(), AtmKeyValStr(), "Numeric"},

Pair this with an implementation type derived from the value variant::

    using GenericAtmKey = decltype(Generic(AtmKeyVal{}))::Const;

The helper gives all metadata users one place to fix the spelling and order
of the atmospheric alternatives. Its returned list must exactly match
``AtmKeyVal``. It currently contains a handwritten string: it centralizes
metadata maintenance but does not automatically track changes to the C++
variant. When adding or reordering a key, update the variant and this helper
together. See the troubleshooting section below if a mismatch reaches the
linker.

Future methods accepting both atmospheric and surface keys can compose their
C++ alternatives with ``Extend``::

    using GenericAtmKey = decltype(Generic(AtmKeyVal{}))::Const;
    using GenericSurfaceKey = decltype(Generic(SurfaceKeyVal{}))::Const;
    using GenericStateKey = GenericAtmKey::Extend<GenericSurfaceKey>;

A corresponding metadata helper must produce the same ordered union:
atmospheric alternatives first, then previously absent surface alternatives.
Do not blindly concatenate strings when the sets overlap: ``Extend`` removes
duplicate types. Both correlation arguments can then use the same combined
key set, allowing either atmospheric/surface ordering. Accepting those types
only establishes the interface; cross-domain covariance mapping and its
validation still need implementation.

Troubleshooting Generic signature and linker errors
--------------------------------------------------

An ``undefined reference`` or ``undefined symbol`` for a workspace method
containing ``Generic<...>`` can mean that its implementation does not match
the generated declaration. C++ treats the mismatching definition as another
overload, so it may compile successfully while leaving the generated method
undefined.

Alternative order is part of the C++ type: ``Generic<Numeric, Vector>`` and
``Generic<Vector, Numeric>`` are different types, even though they accept
the same alternatives. ``Extend`` also preserves order. An alias derived
from a value variant follows that variant's declaration order.

When this happens:

* Find the method declaration in ``build/src/auto_wsm.h`` (adjust the build
  directory if necessary) and compare it with the definition in ``src/m_*.cc``.
* Compare the exact alternative lists and their order with ``gin_type`` and
  ``gout_type`` in the method metadata. Check for missing or extra types,
  pointee constness, and value versus reference parameters as well.
  Top-level ``const`` on a by-value parameter does not change its signature;
  ``const`` on a Generic alternative does.
* Correct the metadata or the implementation, then rebuild to regenerate the
  declaration. Do not edit ``auto_wsm.h`` directly. If the signatures already
  match, check that the source file is included in the build and that the
  definition has the expected namespace and external linkage.

For an earlier compiler diagnostic, enable ``-Wmissing-prototypes`` with
Clang/AppleClang, or ``-Wmissing-declarations`` with GCC. These flags warn
when an externally visible function definition has no matching prior
declaration. Implementation files must include ``workspace.h`` (or the
generated declaration header) for this check to work. With Clang, a Generic
order mismatch then produces a warning at the definition such as::

    warning: no previous prototype for function 'MyMethod' [-Wmissing-prototypes]

Use ``-Werror=missing-prototypes`` with Clang or
``-Werror=missing-declarations`` with GCC to make this a compilation error.
The flags can be added to a single compiler invocation when investigating
an error, or to a developer build's C++ flags. They also diagnose unrelated
externally visible helpers without declarations; give those helpers a
proper declaration or internal linkage as appropriate. Do not follow the
compiler's suggestion to make a workspace method ``static``: the generated
adapter must be able to link to it.
