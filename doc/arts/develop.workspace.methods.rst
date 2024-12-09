Workspace methods
#################

Workspace methods are the core interaction a user has with ARTS calculations.
They are called from the :class:`~pyarts.workspace.Workspace` object in python but are
implemented in C++ under-the-hood.

The workspace methods definition are mainly located in the ``workspace_methods.cpp``
and ``workspace_meta_methods.cpp`` files.  Some workspace methods are 
automatically generated elsewhere.

Workspace methods are exposed to other C++ files via the ``<workspace.h>`` header file.
Only the CMake target ``artsworkspace`` and targets that depend on it may include
``<workspace.h>``.

Defining a workspace method
===========================

You may define either a basic method or a meta-method that depends on other methods.
The easiest way to do so is to copy an existing method and modify it to suit your needs.
Please design your method input and output to allow meta-methods to be defined,
as they simplify composibility signficantly.

Basic methods
-------------

Basic methods completely define their input and output and do not explicitly
depend on other methods.  The method must be manually implemented.

Basic methods are defined entirerly in ``workspace_methods.cpp``.
They are part of a map object called ``wsm_data``.  The name of the
method is the key and the value is a struct with the following fields:

- ``desc`` - a description of the method as a string.
- ``author`` - the author(s) of the method as a list of strings.
- ``out`` - the :doc:`develop.workspace.variables` output of the method as a list of strings.
- ``gout`` - the generic output of the method as a list of strings.  These must not be :doc:`develop.workspace.variables`.
- ``gout_type`` - the type of the generic output as a list of strings.  These must be :doc:`develop.workspace.groups`.
- ``gout_desc`` - the description of the generic output as a list of strings.
- ``in`` - the :doc:`develop.workspace.variables` input of the method as a list of strings.
- ``gin`` - the generic input of the method as a list of strings.  These must not be :doc:`develop.workspace.variables`.
- ``gin_type`` - the type of the generic input as a list of strings.  These must be :doc:`develop.workspace.groups`.
- ``gin_value`` - the default value of the generic input as an optional initialized :doc:`develop.workspace.groups`.
- ``gin_desc`` - the description of the generic input as a list of strings.
- ``pass_workspace`` - a boolean indicating if a :class:`~pyarts.workspace.Workspace` instance should be passed to the method.  If true, the first argument to the method is a ``const Workspace&``.

The expected signature of the method will depend on these fields.
A linker error will likely occur if the actual signature does not match
the expected signature.

The ``in`` and ``out`` may contain the same :doc:`develop.workspace.variables`.  If they do, the variable must be
initialized before the method is called because it is treated as if the method is intended to
simply modify the existing value.  Please indicate strongly in the documentation if you sometimes overwrite the input variable.

On the other hand, if the :doc:`develop.workspace.variables` is only in ``out`` and not in ``in``,
it is treated as if the workspace variable is created by the method.  Note that since the type system
does not account for this, it is important that you clear the current state of the :doc:`develop.workspace.variables`
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

Meta methods
------------

Meta methods do not define all their input and output, but instead define a call
order into other methods.  From this call order, the inputs of the user-facing
workspace method is inferred.  This method should not be implemented manually.

These methods are defined in ``workspace_meta_methods.cpp``.  They are defined
as part of a list called ``wsm_meta``.
A single meta method data contains:

- ``name`` - the name of the method as a string.
- ``desc`` - a description of the method as a string.
- ``author`` - the author(s) of the method as a list of strings.
- ``methods`` - the methods that the meta method depends on as a list of strings.
- ``out`` - the output of the method as a list of strings.  These must be workspace variables.
- ``preset_gin`` - The preset ``gin`` values for the method as a list of workspace values.
- ``preset_gin_value`` - The preset ``gin_value`` values for the method as a list of workspace values.

.. tip::

  A meta method may depend on another meta method.  If it does it is important that the
  meta method it depends on is defined before it in the list.

Automatic methods
-----------------

All methods that execute a workspace agenda are automatically generated.
These will be named as ``agenda_nameExecute`` and may otherwise be
treated as normal workspace method.
You need to do nothing to define these methods.  But please refrain from defining
them manually as that may cause undefined naming conflicts.

The expected signature of the method :func:`~pyarts.workspace.Workspace.propagation_matrix_agendaAuto` is also
generated automatically near the end of ``workspace_methods.cpp``.  It takes
its input and output from a list of other methods.  Feel free to add to this
list but make sure that any naming conflicts regarding ``gin`` are resolved
before doing so.  Adding a method to this list may also require changing the
actual signature (which is why the method is generated, so that a change in
the required actual signature is immediately made apparent).

The methods the begin with ``RetrievalAdd...`` are partly generated.
These methods all require a corresponding ``jacobian_targetsAdd...`` method
that fills in the ``jacobian_targets`` workspace variable.  To keep that
part of the signature consistent, the additional ``RetrievalAdd...`` information
is simply appended to the the ``in``, ``out``, and ``gin``-lists of the
corresponding ``jacobian_targetsAdd...`` method using the local ``jac2ret`` lambda.

Generated files
===============

The workspace method interface generates a lot of files during the build process.
These generated files are located in the build directory and are named
as ``auto_wsm_N.cc``, where N is a number, as ``auto_wsm.cpp``, as ``auto_wsm.h``,
and as ``auto_wsmmeta.cpp`` for the C++ interfacing code.  The python-binding
code is also generated as ``py_auto_wsm_N.cpp``, where N is still a number.

Workspace method naming convention
==================================

Names carry meaning.  Please follow the naming convention below and
please do not hesitate to fix any naming inconsistencies you find.

Method naming
-------------

Workspace method names should be descriptive and follow the naming convention
that the main workspace variable output of the method in ``snake_case``
is followed by a short but descriptive name of what the method does with the output
in ``PascalCase``.
A general rule of thumb is to use verbs for methods that modify the workspace
variable and nouns for methods that create a new workspace variable.

For example, :func:`~pyarts.workspace.Workspace.propagation_matrixAddLines`
has a main output of :attr:`~pyarts.workspace.Workspace.propagation_matrix` and
adds line absorption to it.  It needs to be preceeded by a call to 
:func:`~pyarts.workspace.Workspace.propagation_matrixInit` which sets up the
propagation matrix to an initial state.

Of course, every use-case is different but please try to follow this convention.

File naming
-----------

The file that a workspace method is implemented in should be named ``m_<concept>.cc``.
The concept should be a short but descriptive name of what the methods therein do.
Multiple methods per file is allowed and encouraged, but keep them conceptually similar.
To ensure compatibility with various filesystems, please avoid using spaces
and capital letters in the filename.

Lastly, please ensure that the file is listed in the CMake target ``artsworkspace``
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
