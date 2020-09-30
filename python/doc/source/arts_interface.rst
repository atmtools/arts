The interactive ARTS interface
==============================

The interactive ARTS interface is implemented by the
:class:`pyarts.workspace.Workspace` class. Each object of this
class represents an ARTS workspace on which ARTS workspace methods
can be executed:

.. code-block:: python

    from pyarts.workspace import Workspace
    ws = Workspace()

Each Workspace object exposes all available ARTS workspace variables
(WSVs) and workspace methods (WSMs) as attributes:

.. code-block:: python

   >>> ws.y
   ARTS Workspace Variable: y
   >>> ws.yCalc
   ARTS Workspace Method: yCalc

Documentation for WSVs or WSMs can be printed by calling their
:code:`help()` member functions:

.. code-block:: python

   >>> ws.y.help()
   ARTS Workspace Variable

   Name:  y
   Group: Vector

   The measurement vector.
   [...]

Calling workspace methods
-------------------------

To execute a WSM on a workspace method it suffices to call the
corresponding member function of a workspace instance. For example
to call the :code:`yCalc` WSM:

.. code-block:: python

   ws.yCalc()

The function-call syntax is the same as in an ARTS controlfile. There
are two ways to provide function arguments to a WSM call:

1. Using positional arguments:

.. code-block:: python

   ws.yCalc(ws.yf)


2. Using named arguments:

.. code-block:: python

   ws.yCalc(y = ws.yf)

Both calls replace the output :code:`y` of the :code:`yCalc` workspace methods
with the :code:`yf` WSV. Similar as in a controlfile, arguments that are not
listed are replaced with their defaults. 

   
Accessing workspace variables
-----------------------------

The WSV attributes of a :class:`Workspace` object provide a symbolic
representation of the variable. Therefore they can only have a value when they
are associated to a specific workspace.

Setting the value of a workspace variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The value of a workspace variable in a given workspace can be set by
assigning directly to the corresponding attribute of the :class:`Workspace`
object:

.. code-block:: python

   >>> ws.y = np.ones(10)

When a value inside an ARTS workspace is set, the incoming data is always
copied. This means that when an Python array variable is assigned to a workspace
variable, the array variable can be in-place modified without changing the value
of the ARTS WSV.

.. code-block:: python

   >>> x = np.ones(10)
   >>> ws.y = x
   >>> x[0] = 0.0
   >>> x[0] == ws.y.value[0]
   False

To assign a value to an ARTS WSV, it has to be compatible with
the corresponding ARTS group. The table below summarizes the
mapping of ARTS groups to Python types.

   +----------------------------------+---------------------+
   | ARTS group                       +  Python type        |
   +----------------------------------+---------------------+
   | Index                            +  :code:`int`        |
   +----------------------------------+---------------------+
   | Numeric                          +  :code:`float`      |
   +----------------------------------+---------------------+
   | Vector, Matrix, Tensor[3,...,7]  +  :code:`numpy.array`|
   +----------------------------------+---------------------+
   | ArrayOfIndex                     +  :code:`list`       |
   +----------------------------------+---------------------+
   | ArrayOfString                    +  :code:`list`       |
   +----------------------------------+---------------------+
   | Sparse                           + :code:`scipy.sparse`|
   +----------------------------------+---------------------+

In addition to the groups, the **pyarts** package provides a number of
specialized classes to represent ARTS groups. Refer to :ref:`ARTS classes` for
an overview.

.. note:: The interface performs some simple conversions in order to simplify
   assigning values to WSVs. For groups :code:`Vector`, :code:`Matrix` and
   :code:`Tensor`, the assigned value is casted to the corresponding shape,
   which may lead to unexpected results when for example setting the
   :code:`z_field` with a 1D :code:`numpy.array`.

Accessing the value of a workspace variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The value of a workspace variable within a given workspace can be accessed
through its :code:`value` attribute:

.. code-block:: python

   >>> ws.y.value
   array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])

WSVs belonging to any of the groups :code:`Vector`, :code:`Matrix`,
:code:`Tensor[3, ..., 7]` are passed back from ARTS as references. Their
values can therefore be changed by manipulating the returned object:

.. code-block:: python

   >>> y = ws.y.value
   >>> y[0] = 0.0
   >>> ws.y.value
   array([0., 1., 1., 1., 1., 1., 1., 1., 1., 1.])

.. note:: Since all array-type variables which hold numeric data are
   returned from ARTS by reference, extracting results from a simulation
   should be done by copying the variable using e.g. :code:`numpy.copy`.
   Otherwise the values will be overwritten when new calculations are
   performed on the workspace.

Executing control files
-----------------------

:code:`*.arts` controlfiles can be executed on a workspace using the
:code:`execute_controlfile` member function of the workspace object. For
example, to prepare a workspace using the utility controlfiles distributed
with ARTS:

.. code-block:: python

   ws.execute_controlfile("general/general.arts")
   ws.execute_controlfile("general/continua.arts")
   ws.execute_controlfile("general/agendas.arts")
   ws.execute_controlfile("general/planet_earth.arts")

Defining and executing agendas
------------------------------

It is also possible to specify ARTS agendas in Python. This is done by defining
a suitable function and using the :code:`@arts_agenda` decorator to transform it
into an agenda.

.. code-block:: python

   from pyarts.workspace import arts_agenda

   @arts_agenda
   def ppath_agenda(ws):
       ws.Ignore(ws.rte_pos2)
       ws.ppathStepByStep()

   # Copy ppath_agenda into workspace.
   ws.ppath_agenda = ppath_agenda

The format required for turning a Python function into an ARTS agenda is that it
takes a single input argument, :code:`ws` in the example, which represents the
workspace. ARTS WSMs can be executed in the agenda calling them on the functions
input argument.

INCLUDE statements
^^^^^^^^^^^^^^^^^^

For consistency with ARTS controlfile syntax, also :code:`INCLUDE` directives
are supported by the Python interface. The :code:`INCLUDE` directive is used to
include a controlfile (or another agenda defined in Python) in an agenda
definition.

.. code-block:: python

   from pyarts.workspace import arts_agenda

   @arts_agenda
   def ppath_agenda(ws):
       INCLUDE("my_controlfile.arts")


Python within ARTS
^^^^^^^^^^^^^^^^^^

It is even possible to execute Python code within an agenda. By replacing the
ARTS :code:`iy_space_agenda`, this could for example be used to perform
simulations which assume a different cosmic microwave background temperature:

.. code-block:: python

        import scipy.constants as c

        @arts_agenda
        def space_agenda(ws):
            # Since everything happens in Python we need
            # to tell ARTS that we are using all in and outputs.
            ws.Ignore(ws.f_grid)
            ws.Ignore(ws.rtp_pos)
            ws.Ignore(ws.rtp_los)
            ws.Touch(ws.iy)

            # Temperatures and frequency
            t = 4.735 # Some men just want to watch the world burn.
            f = ws.f_grid.value

            # Compute radiances
            c1 = 2.0 * c.h / c.c ** 2
            c2 = c.h / c.k
            b = c1 * f ** 3 / (np.exp(c2 * f / t) - 1.0)

            # Put into iy vector.
            ws.iy = np.zeros((f.size, ws.stokes_dim.value))
            ws.iy.value[:, 0] = b

        # Copy ppath_agenda into workspace.
        ws.iy_space_agenda = space_agenda
