Workspace options
#################

Workspace options are wrappers for creating ``enum class`` types that fully define the set of valid options 
that are required for "some" task.  They are used instead of strings to ensure that the methods that use them
always cover the full set of solutions required (i.e., by switching over all options in a switch-statement).

They are mostly intended as input to workspace methods, but they may also be part of variuous workspace groups
to define the logic of the group.  The options are defined in the ``arts_options.cc`` file as parts of the
``opts`` list object.  Each option is a struct with the following fields:

- ``name`` - the name of the option as a string.
- ``desc`` - a description of the option as a string.
- ``values_and_desc`` - a list of lists of strings.  The outer-most list defines a unique option.  The inner lists must all have the same length.  The inner list defines at least two strings.  The first value of the inner list are used to build the ``enum class`` value and must thus be usable as valid C++ code.  The last value of the inner list is the description of the option.  All values inbetween are various ways to spell the option.  The ``enum class`` will be constructible from strings of any of these values.

What do you get by defining an option?
======================================

Using ``T`` as the option type, ``e`` as an instance of the option type, and ``s`` as an instance of ``std::string``,
you will automatic have access to the following interface for each option type:

- ``good_enum(e)`` - returns a boolean that checks if the option is valid.
- ``toString(e)`` - returns a string representation of the option.  Note that this is a templated type and you can input a templated index to select the name you want from the lists of strings.
- ``to<T>(s)`` - returns the option from a string.  It will look through all provided names for the option type.
- ``operator<<(std::ostream &os, const T x)`` - allows streamed output of the option.
- ``operator>>(std::istream &is, T &x)`` - allows streamed input of the option.
- In addition to this, all the requirements for being a :doc:`develop.workspace.group` are also automatically fulfilled.

Generated files
===============

The workspace group interface generates ``enums-common-helper.h``, ``enums.cpp``, ``enums.h``, and ``enumsNAME.h`` files for the C++ interface, here ``NAME`` is the name of the option.
The Python binding is generated as the ``py_auto_options.cpp`` and ``py_auto_options.h`` files.

Workspace variable naming convention
====================================

Workspace options may be named as pleased.  The name should be descriptive of the type of options that it offers.
