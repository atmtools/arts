The Workspace
=============

The workspace is the primary interface to interact with the pyarts
module. It is a collection of rules relating groups, variables, methods, and agendas.

Workspace groups are the classes that can operate on the workspace.  These are the
logic behind what type of data is passed around in the model.

Workspace variables are the named instances of these workspace groups.
There are two types of workspace variables: named and user-defined.
Named variables are defined by the workspace and are always considered available.
A variable is initialized once written to.
This can be done either by assigning a
value to the workspace variable, or by calling a workspace method that gives that variable as an output.

Workspace methods are functions that perform some task.
The input and output of these methods are either workspace variables or generic instances of workspace groups.
If a named workspace variable has been initializes on the workspace, it is not necessary
to specify it in the interface when calling a method requiring that variable.
As a simple example, if ``ws.A()`` takes ``a`` as a named input, and ``ws.a`` has been initialized,
then ``ws.A()`` can be called without specifying ``a`` as an input.  The implicit call here is 
``ws.A(a=ws.a)``.  It is possible to also define input generically, like ``ws.A(a=1)``,
assuming the workspace group of ``a`` can be initialized from the value ``1``.
Unlike python in general, methods have pure output variables as is more common
in ``C++``.  Note that this means that methods generally modify the state of the
workspace when called.

Agendas are collections of methods that are executed to perform some task inside
other methods.  Think of them as callback methods that are executed within the
context of other methods.  Agendas always have access to all named variables at
the level they are executed at, but when executed inside a method agendas will not
affect the state of workspace variables outside its own scope. 
Agendas make it possible to define very complicated
methods and concepts in a modular way.  See the Examples in the documentation
for more information.

.. toctree::
    :hidden:

    workspace.groups
    workspace.methods
    workspace.variables
    workspace.agendas

.. rubric:: Workspace Groups
.. include:: workspace.groups.auto.rst

.. rubric:: Workspace Variables
.. include:: workspace.variables.auto.rst

.. rubric:: Workspace Methods
.. include:: workspace.methods.auto.rst

.. rubric:: Workspace Agendas
.. include:: workspace.agendas.auto.rst
