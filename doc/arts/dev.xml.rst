XML IO
======

This is the main way that we perform file IO in ARTS.

Each workspace group must support XML file IO.
If a group is added that does not support XML file IO, ARTS cannot compile.
Support is done by overloading a templated struct, ``xml_io_stream``,
that defines the core concepts involved in the workflow.

How the group handles XML file IO inside this is up to the group itself.
Generally, all lower level overloads should already be defined, so it
is common to just chain several calls to other XML IO methods.

The ``xml_io_stream``
---------------------

This is the core struct as defined in ``xml_io_struct.h``

.. code-block:: C++


