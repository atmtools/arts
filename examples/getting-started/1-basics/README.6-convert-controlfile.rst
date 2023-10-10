To take a legacy ARTS controlfile (``.arts`` file) as a starting point, you can use the ``arts_convert.py`` script to convert it to a Python script. The script is included with the conda package and should be available in your path. For example:

.. code-block:: bash

    arts_convert.py mycontrolfile.arts

This will create a Python script called ``mycontrolfile.py`` which you can then edit and run.

.. warning::
    Before converting your controlfile, make sure it is compatible with the current version of ARTS. Otherwise, the conversion may fail.

ARTS Controlfile:

.. code-block:: ruby

    Arts2 {
        VectorNLinSpace( f_grid, 5, 320e9, 322e9 )
    }

Python conversion:
