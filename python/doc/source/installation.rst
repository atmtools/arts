Installation
============

Pre-compiled binaries
^^^^^^^^^^^^^^^^^^^^^

Pre-compiled binaries of the ARTS interface for macOS and Linux including the
ARTS engine for performing RT simulations can be installed in a `Miniforge3
<https://github.com/conda-forge/miniforge#miniforge>`_ environment from the
rttools channel:

.. code-block:: bash

    mamba install -c rttools pyarts

Other versions are available on the `ARTS homepage <https://radiativetransfer.org/getarts/>`_.

ARTS depends heavily on catalog data. It is recommended to call :py:func:`pyarts.cat.download.retrieve` at the beginning of your Python scripts to download and cache the latest version of the ``arts-cat-data`` and ``arts-xml-data`` packages. Alternatively, you can download the matching catalog data manually from the `Github release <https://github.com/atmtools/arts/releases/>`_ page.

.. warning::
    The pyarts package was created for use with Miniforge3.  If you are using
    Anaconda, you will need to install it into a separate environment with the
    ``conda-forge`` channel:

    ``conda create -n arts -c conda-forge -c rttools pyarts``

    But compatiblity issues remain, e.g. it has been reported that the package
    fails to install on Intel Macs with Anaconda. Therefore, we recommend to
    use Miniforge3.


Building from source
^^^^^^^^^^^^^^^^^^^^

See top level ``README.md`` at `our github repository <https://github.com/atmtools/arts>`_.

The API of ARTS is constantly being updated to suit our and our users future need.
Do not expect code compatibility between any two patches.

We strongly recommend to use the latest version of ARTS from conda unless you want to help with the development.
