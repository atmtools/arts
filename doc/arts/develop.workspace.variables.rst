Workspace variables
###################

Workspace variables allow named data to be passed into workspace methods.
They are attributes of the :class:`~pyarts.workspace.Workspace` object.

Most workspace variables are defined in the ``workspace_variables.cpp`` file.
Workspace agendas are automatically defined as workspace variables from the ``workspace_agendas.cpp``.

Defining a workspace variable
=============================

In ``workspace_variables.cpp``
------------------------------

The workspace variable definitions are located in the ``workspace_variables.cpp``
as part of the ``wsv_data`` map object.  The name of the variable is the key
of the map and the object are a struct with the following fields:

- ``desc`` - a description of the variable as a string.
- ``type`` - the workspace group of the variable as a string.
- ``default_value`` - an optional default value of the workspace variable that is defined when a :class:`~pyarts.workspace.Workspace` object is default-initialized.

In ``workspace_agendas.cpp``
----------------------------

See the page on :doc:`develop.workspace.agendas` for more information.  The workspace variable defined from the workspace agenda will
have no default value.  It will be of the type :class:`~pyarts.arts.Agenda` by default but of type :class:`~pyarts.arts.ArrayOfAgenda`
if the field ``array`` is ``true``.

What qualifies as a workspace variable?
=======================================

The intent of a workspace variable is to pass information between workspace methods.
Workspace variables must therefore be both input and output of different workspace methods.

If a method has output that is never used by a different method, it should not be a workspace variable. Instead, it should be a local variable of the method.  As defined by workspace method, it should be a ``gout``.

If a method has input that is not produced by a different method, it should not be a workspace variable. Instead, it should be a local variable of the method.  As defined by workspace method, it should be a ``gin``.

There are only a few exceptions to this rule:

- Workspace agendas are not required to be output of a workspace method.  They are expected to be consumed by a workspace method.  The reason for this exception is that they define a special logic that may not be easily expressed as a workspace method or a workspace option.
- Workspace agendas' inputs and outputs are also not required to follow the above rule.  The reason is the same as for why workspace agendas are themselves an exception.
- Workspace options are not required to be output of a workspace method.  The reason is that workspace methods must fully define their logic, and workspace options are the explicitly named options that a workspace method is expected to understand.

Generated files
===============

The workspace variable interface generates ``auto_wsv.cpp`` and ``auto_wsv.h``
for the C++ interface to the workspace variables.  The Python binding is generated
as the ``py_auto_wsv_N.cpp`` files, where N is an integer.

Workspace variable naming convention
====================================

Names carry meaning.  Please follow the naming convention below, and
please do not hesitate to fix any naming inconsistencies you find.

Workspace variables should be named in ``snake_case``.  The name should be
explicit.  Avoid TLAs [#f1]_ and other abbreviations that are not world-wide exclusive.
A common TLA to avoid is ``lte`` because it means different things in different
contexts (e.g., local thermodynamic equilibrium or long-term evolution).
Using ``lte`` in a name is therefore ambiguous and strongly discouraged.  On-the-other hand, ``nonlte``
is unambiguous - it is not used as anything other than non-local thermodynamic equilibrium.

To allow searching and connecting different workspace variables to their intended context,
we have reserved a set of terms that should be used in the names of related workspace variables.
This list will connect workspace variables under a rubric ``Related workspace variables`` in the documentation.
Please note that combining more than one of these terms is possible and highly encouraged.
The workspace variable :attr:`~pyarts.workspace.Workspace.propagation_matrix` along the :attr:`~pyarts.workspace.Workspace.ray_path`
is called :attr:`~pyarts.workspace.Workspace.ray_path_propagation_matrix` to connect it to both concepts.
The full list is available in the ``workspace_variables_keywords.cpp`` file.
As of writing this, the list includes the following terms:

- ``absorption`` - for workspace variables related to absorption.
- ``jacobian`` - for workspace variables related to Jacobian.
- ``spectral_radiance`` - for workspace variables related to spectral radiance.
- ``propagation_matrix`` - for workspace variables related to propagation matrices.
- ``source_vector`` - for workspace variables related to source vectors.
- ``nonlte`` - for workspace variables related to non-local thermodynamic equilibrium.
- ``ray_path`` - for workspace variables related to ray paths.
- ``scattering`` - for workspace variables related to scattering.
- ``grid`` - for workspace variables that are also grids.
- ``measurement`` - for workspace variables related to measurements.
- ``model_state`` - for workspace variables related to model states.
- ``disort`` - for workspace variables related to DISORT.

Check the ``workspace_variables_keywords.cpp`` file for the correct list regardless.  If you read this and the list above is no longer correct, please update it.

.. [#f1] TLA means three-letter abbreviations, and it is commonly used in text and speech to indicate why you should avoid using TLAs - abbreviations feel easy to use and common-know-how once you get used to them but the first time you encounter one you have no idea what it means, which distracts from the content being discussed.
