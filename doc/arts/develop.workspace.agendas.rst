Workspace agendas
#################

At their simplest, workspace agendas uses a set of :doc:`develop.workspace.methods` that combined requires
a set of :doc:`develop.workspace.variables` to produce another set of workspace variables.

Workspace agendas are defined in the ``workspace_agendas.cpp`` file.
The workspace agenda definitions are located in the ``wsa_data`` map object.  The name of the agenda is the key of the map and the object is a struct with the following fields:

- ``desc`` - a description of the agenda as a string.
- ``outputs`` - the workspace variables that are produced by the agenda as a list of strings.
- ``inputs`` - the workspace variables that are consumed by the agenda as a list of strings.
- ``enum_options`` - a list of strings that are the options for the agenda.  Defaults to empty.
- ``enum_defaults`` - the default enum option for the agenda.  Defaults to no-default by being empty.
- ``array`` - a boolean for whether the agenda is described as an :class:`~pyarts.arts.ArrayOfAgenda` (if true) or an :class:`~pyarts.arts.Agenda` (if false, default).

The two enum-variables, if set, will generate workspace methods that allow setting the agenda to a specific named option.
This is useful for agendas that have known "modes".  However, it puts some constraints on the implementation of the agenda.
The implementor must implement a method ``get_MYAGENDA(const std::string&)`` in the
``workspace_agenda_creator.h`` and ``workspace_agenda_creator.cpp`` files.
Please see existing options in those files for examples.

What qualifies a workspace agenda?
==================================

A method or a set of methods must be able to use the workspace agenda to produce a set of desired workspace variables.
Generally, if you can keep the number of inputs limited and the number of outputs fixed, you may consider a workspace agenda.
However, overpopulating the workspace with agendas is not recommended. So please reuse existing ones as appropriate.

Take :attr:`~pyarts.workspace.Workspace.ray_path_observer_agenda` as an example.  The different ways we can define a path
from an observer to a background is practically limitless.  However, if we are interested in :attr:`~pyarts.workspace.Workspace.spectral_radiance`
for some observer and frequency,
we know that we just need the background radiation from the background of the observed :attr:`~pyarts.workspace.Workspace.ray_path` and the state of the 
atmosphere along the path.  So a method like :func:`~pyarts.workspace.Workspace.measurement_vectorFromOperatorPath`, which takes a set of frequencies and a set of observers,
may effectively ignore the complexity of computing the :attr:`~pyarts.workspace.Workspace.ray_path` internally, and leave that to the workspace agenda.

Generated files
===============

The workspace group interface generates ``auto_wsa.cpp`` and ``auto_wsa.h`` files for the C++ interface.
No additional files are generated for the Python binding.

Workspace agenda naming convention
==================================

Workspace agendas should be named in ``snake_case`` as they are also :doc:`develop.workspace.groups`.

Do not hesitate to implement helpers
====================================

Several agendas have helper methods to make the user interface easier.
These methods should follow a simple naming convention.

- There are methods that set an agenda to a named option, such as :func:`~pyarts.workspace.Workspace.spectral_radiance_space_agendaSet`.
- There are methods that compute an agenda based on existing :doc:`develop.workspace.variables`, such as :func:`~pyarts.workspace.Workspace.propagation_matrix_agendaAuto`.
- There are methods that set and agenda based on multiple options, such as :func:`~pyarts.workspace.Workspace.Setup`.
