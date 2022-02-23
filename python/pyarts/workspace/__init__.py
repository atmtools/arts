"""

The `workspace` subpackage provides an interactive interface to ARTS and
can be used to directly call workspace methods and agendas as well as
access and manipulate workspace variables.

Setup
-----

The CMake build process should automatically generate this interface as long as
you build the pyarts package.  A key here is that an internal representation of
the builtin Arts functions and classes are generated using the pybind11 library
for the bindings.

The Workspace Class
-------------------

The main functionality of the interface is implemented by the Workspace class.
A Workspace object represents an ongoing ARTS simulation and is used to execute
controlfiles and workspace methods and access workspace variables

>>> from arts.workspace import Workspace
>>> ws = Workspace()

For a basic interface, see

>>> help(ws)

Executing Controlfiles
----------------------

Controlfiles can be executed on a workspace using the
:func:`Include` method.

>>> Include(ws, "general/general.arts")

The search path for controlfiles is the current directory, the paths
provided in the environment variable ``ARTS_INCLUDE_PATH``, and the search paths
that were provided in the ``ARTS_DEFAULT_INCLUDE_DIR`` variable while building
Arts.  Controlfiles are parsed inline, so repeated parsing of a controlfile
is possible.  If this fails, the controlfile is effectively not parsed.  (A key
exception here is that variables that are created by a controlfile are
effectively done so by the parsing and not during the execution of the
controlfile.  As it is not allowed to create a named variable twice in a
controlfile, this means that a controlfile that creates variables often cannot
be included twice, regardless of intitial success.)

Calling Workspace Methods
-------------------------

ARTS workspace methods are available as member functions of each Workspace
object:

>>> ws.AtmosphereSet1D()
>>> ws.IndexSet(ws.stokes_dim, 1)

Arguments can be passed to a workspace function in three ways:

    1. As workspace variables using the attributes of the
       workspace object, such as :code:`ws.stokes_dim` in the
       example above.
    2. Using one of classes available in
       :mod:`pyarts.classes`
    3. Passing supported python objects directly (limited to pure input)

Arguments to a WSM can be passed using either positional or named arguments or
both.

>>> ws.VectorNLogSpace(ws.p_grid, 361, 500e2, 0.1 )

Available keywords are the names of the generic outputs and inputs defined in
methods.cc

>>> ws.abs_speciesSet(species=[ "O3" ])

The documentation is generally available via python's help() command, e.g.,

>>> help(ws.yCalc)

Calls to supergeneric functions are resolved by the interface.

Workspace Variables
-------------------

Variable objects are associated to a workspace, and can be accessed as
attributes of a workspace, i.e. using for example:

>>> ws.y

their value can be accessed using the value property as

>>> ws.y.value

This will return a builtin representation of the value of the workspace if it is
initialized.  You should be able to set workspace variables from most python ways
of representing the type such as

>>> ws.y = [1, 2]

or

>>> ws.y = np.array([1,2])

This will initialize the workspace variable.

For more information, see for example

>>> help(ws.y)

or

>> help(ws.y.value)

"""

from pyarts.workspace.workspace import Workspace, arts_agenda, Include # noqa 
