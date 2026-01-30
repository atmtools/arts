Installation
============

Pre-compiled binaries
^^^^^^^^^^^^^^^^^^^^^

Pre-compiled binaries of the ARTS interface for macOS and Linux including the
ARTS engine for performing RT simulations can be installed in a `Miniforge3
<https://github.com/conda-forge/miniforge#miniforge>`_ environment from the
rttools channel:

.. code-block:: bash

    mamba install -c rttools-dev pyarts3

.. warning::
    The pyarts package was created for use with Miniforge3.  If you are using
    Anaconda, you will need to install it into a separate environment with the
    ``conda-forge`` channel:

    ``conda create -n arts -c conda-forge -c rttools-dev pyarts3``

    But compatiblity issues remain, e.g. it has been reported that the package
    fails to install on Intel Macs with Anaconda. Therefore, we recommend to
    use Miniforge3.


Building from source
^^^^^^^^^^^^^^^^^^^^

See top level ``README.md`` at `our github repository <https://github.com/atmtools/arts>`_.
By default, this shows ARTS 3.  Use the dropdown menu to
select branch or tagged versions if you want to build ARTS 2.6 or older.

ARTS 3 is under active development.
Do not expect code compatibility between any two patches.
