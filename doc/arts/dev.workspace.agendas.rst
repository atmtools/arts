Workspace agendas
#################

At their simplest, workspace agendas uses a set of :doc:`dev.workspace.methods` that combined requires
a set of :doc:`dev.workspace.variables` to produce another set of workspace variables.
All agendas have access to the global workspace but are expected to only directly require
a limited subset to produce their outputs.

Workspace agendas are defined in the ``workspace_agendas.cpp`` file.
The workspace agenda definitions are located in the ``wsa_data`` map object.
The name of the agenda is the key of the map.  Throughout this text, we will use ``agendax``
to represent a generic agenda name.
The agenda object itself is a struct with the following fields:

- ``desc`` - a description of the agenda as a string.
- ``outputs`` - the workspace variables that are produced/output by the agenda as a list of strings.
- ``inputs`` - the workspace variables that are consumed by the agenda as a list of strings.
- ``enum_options`` - (optional) a list of strings that are the options for the agenda.  Defaults to empty.
- ``enum_defaults`` - (optional) the default enum option for the agenda.  Defaults to no-default by being empty.
- ``output_constraints`` - (optional) a vector of structs that define constraints on the output after the agenda has been executed.  Defaults to unconstrained.

The two enum-variables, if set, will generate workspace methods that allow setting the agenda to a specific named option.
This is useful for agendas that have known "modes".  However, it puts some constraints on the implementation of the agenda.
The implementor must implement a method ``get_agendax(const std::string&)`` in the
``workspace_agenda_creator.h`` and ``workspace_agenda_creator.cpp`` files.
Please see existing options in those files for examples.

If you want to ensure that the outputs have the right size or correct values (e.g., positive values),
after executing the agenda, you can use the ``output_constraints`` field.  This is a complicated field
that is a vector of constraints.  Each constraint has a test and a message.  The test should be 
a compilable expression that returns a boolean true when the constraint is satisfied.  The message
is the string that is printed in the exception when 1) the test fails, and 2) on the documentation page.
Additionally, any number of expressions can be added to the message.  These expressions added to the message
and evaluated at runtime when a failure occurs (e.g., ``{x.size()}``).

How to add a workspace agenda?
==============================

The easiest way to add an agenda is to copy an existing agenda.
Go to the ``workspace_agendas.cpp`` file.
Copy an existing agenda and modify it to fit your needs.

If you want to provide user-friendly settings for your agenda, use the ``enum_options`` and ``enum_defaults`` fields.
You will need to add the ``get_agendax(const std::string&)`` method in the ``workspace_agenda_creator.h`` and
``workspace_agenda_creator.cpp`` files if you do so, or you will encounter a compile time error.
Adding these methods is once again easy to do by copying an existing method and modifying it to fit your needs.

.. note::

  Keep the switch-statement structure of the ``get_agendax(const std::string&)`` methods.
  All modern compilers will warn when options are missing from the switch statements.
  Do not solve this warning by adding a ``default`` case to the switch statement as this
  will mute all future warnings.  All cases should be explicitly handled.

What qualifies a workspace agenda?
==================================

A method or a set of methods must be able to use the workspace agenda to produce a set of desired workspace variables.
Generally, if you can keep the number of inputs limited and the number of outputs fixed, you may consider a workspace agenda.
However, overpopulating the workspace with agendas is not recommended. So please reuse existing ones as appropriate.

Take :attr:`~pyarts3.workspace.Workspace.ray_path_observer_agenda` as an example.  The different ways we can define a path
from an observer to a background is practically limitless.
However, if we are just interested in :attr:`~pyarts3.workspace.Workspace.spectral_radiance`
for some observer and frequency,
we know that we just need the background radiation from the background of
the observed :attr:`~pyarts3.workspace.Workspace.ray_path` and the state of the 
atmosphere along the path.
So a method like :func:`~pyarts3.workspace.Workspace.measurement_vectorFromOperatorPath`,
which takes a set of frequencies and a set of observers,
may effectively ignore the complexity of computing
the :attr:`~pyarts3.workspace.Workspace.ray_path` internally, and leave that to the workspace agenda.

Generated files
===============

There are several files generated for both the pure C++ and for the python interface.

Check these if anything is unclear after reading this documentation to
ensure there are no bugs in the generated code.

These are the pure C++ files:

- ``auto_wsa.h``: Contains recursive agenda call handling logic. Also defines the workspace methods ``agendaxExecute`` and ``agendaxSet`` interfaces.
- ``auto_wsa.cpp``: Creates a workspace variable with the name of the agenda.  Implements what is in the namesake header file.  Also implements the logic for documenting these methods.
- ``auto_wsa_options.h``:  Implements the options used in set-methods.  These are stripped down versions of the standard ARTS options.
- ``auto_agenda_operators.h``:  Defines a callable workspace group to represent the agenda.  The type is called ``agendaxOperator``.
- ``auto_agenda_operators.cpp``:  Implements workspace methods ``agendaxExecuteOperator`` and ``agendaxSetOperator``.  The former executes the workspace group from the header file, the latter allows setting the agenda to call a custom operator, e.g., a python function.

This is the only python interface files:

- ``py_auto_agenda_operators.cpp``:  Implements the GIL-safe ``agendaxOperator`` callback types.  It should be possible to create these using any python callable, such as a ``def`` method or a class that implements ``__call__``.  The error should be propagated as if it were in python.

Workspace agenda naming convention
==================================

Workspace agendas should be named in ``snake_case`` as they are also :doc:`dev.workspace.groups`.

Do not hesitate to implement helpers
====================================

Several agendas have helper methods to make the user interface easier.
These methods should follow a simple naming convention.

- There are methods that set an agenda to a named option, such as :func:`~pyarts3.workspace.Workspace.spectral_radiance_space_agendaSet`.  These are default generated by using the enum-variables.
- There are methods that compute an agenda based on existing :doc:`dev.workspace.variables`, such as :func:`~pyarts3.workspace.Workspace.spectral_propmat_agendaAuto`.
- There are methods that set an agenda based on multiple options, such as :func:`~pyarts3.workspace.Workspace.disort_settings_agendaSetup`.
