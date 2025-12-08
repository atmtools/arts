Workspace groups
################

Workspace groups are the types of the :doc:`dev.workspace.variables`.
They define the logic that can be applied to the data.
Workspace groups are all available in the :py:mod:`~pyarts3.arts` module.

Most workspace groups are defined in the ``workspace_groups.cpp`` file.
Additionally, :doc:`dev.workspace.options` are automatically turned into workspace groups using their definitions in the ``arts_options.cpp`` file.

How to add a workspace group?
=============================

First create and name the C++ struct, class, or type define in a header file.
After that, the easiest way to create a new workspace group is to simply copy
an existing workspace group in ``workspace_groups.cpp`` and modify it as
necessary.  Ensure that the compilation of the ``artsworkspace`` target can
see the header file that contains the definition of the workspace group.

You will also need to create the python bindings for the new workspace group.
Again, the easiest way to do this is to simply copy an existing python interface.
It is then possible to overload as many methods and operators as required to make
your new workspace group interface smooth from python.

Note that there are several templated methods available that will help your
class have a better python interface.  Please write more if you find some
defaults are missing.  The most important one for all workspace groups is
``generic_interface``, which handles bindings for all the standard
interfaces required for workspace groups (such as IO).

Defining a workspace group
==========================

In ``workspace_groups.cpp``
---------------------------

The workspace group definitions are located in the ``workspace_groups.cpp``
as part of the ``wsg_data`` map object.  The name of the group is the key
of the map and the object is a struct with the following fields:

- ``file`` - the main header file that must be included to use the workspace group.
- ``desc`` - a description of the variable as a string.
- ``value_type`` - a boolean that defines whether a workspace group instance of this type is copyable in python.  You cannot copy python types such as ``int`` and ``str``.
- ``map_type`` - a boolean that defines whether the workspace group is a map.  Maps can be accessed via the ``[]`` operator using the key type.  Remember that both the ``mapped_type`` and the ``key_type`` of the map should also be workspace groups.

What qualifies as a workspace group?
====================================

You need to ensure that the following code is possible for each workspace group ``T``.
Read this assuming that ``a`` is an instance of ``T``, ``x`` is an instance of
``std::formatter<T>``, ``i`` is an instance of an integer, and ``k`` is an instance of a "key".
The following must compile:

- ``xml_read_from_stream(std::istream &, T&, bifstream *)`` - allow reading the workspace group from an XML file.  Failure to comply leads to a linker error.
- ``xml_write_to_stream(std::ostream &, const T&, bofstream *, const String &)`` - allow writing the workspace group to an XML file. Failure to comply leads to a linker error.
- ``T{}`` - allow default construction.  Failure to comply leads to a static assertion as the group fails the ``WorkspaceGroupIsDefaultConstructible`` concept.
- ``T{a}`` - allow copy construction of.  Failure to comply leads to a static assertion as the group fails the ``WorkspaceGroupIsCopyable`` concept.
- ``std::format("{}", a)``, ``std::format("{:sqNB,}", a)``, ``x.inner_fmt().tags``, and ``static_assert(std::same_as<format_tags, std::remove_cvref_t<decltype(x.inner_fmt().tags)>>)`` - allow formatting the group to a string.  The exception to this are classes that pass one of these consepts: ``std::integral<T>`` or ``std::floating_point<T>`` or ``std::same_as<T, std::string>``.   Failure to comply leads to a static assertion as the group fails the ``arts_formattable_or_value_type`` concept.  Note that to ensure this, you should read the :doc:`dev.classes.formatter` documentation.
- ``x[i]`` and ``x[k]`` should also implement all the above if the group is an array or map type, respectively.  This also holds true for the group of ``k`` for map types.  Failure to compile will lead to difficult errors in the python binding compilation.

There exist a helper method ``generic_interface`` you can send a mutable version
of the ``py::class_<T>`` object to that should ensure that you pass most of the standard tests.

Generated files
===============

The workspace group interface generates ``auto_wsg_init.cpp``, ``auto_wsg_share.h``,
``auto_wsg.h`` and ``auto_workspace.cpp`` files for the C++ interface.
It also generates ``py_auto_wsg.cpp``, ``py_auto_wsg.h``, and ``py_auto_wsgdocs.h`` for the Python binding.

Workspace group naming convention
=================================

Workspace groups should be named in ``PascalCase``.  The name should be
descriptive of the logic that is handled by the group.

File naming
===========

This is up to the developer.  However, try to place the group as early in the compilation as possible.
Preferably, the group should only be part of the ``artscore`` CMake target (or one of its dependencies).
