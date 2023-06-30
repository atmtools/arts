Installation
============

Pre-compiled binaries
^^^^^^^^^^^^^^^^^^^^^

Pre-compiled binaries of the ARTS interface including the ARTS engine for
performing RT simulations can be installed from the conda-forge environment
at the rttools channel:

.. code-block:: bash

    conda install -c rttools pyarts

Building from source
^^^^^^^^^^^^^^^^^^^^

See top level ``README.md`` at `our github repository <https://github.com/atmtools/arts>`_.

Older versions of Arts are available as well to support the reproduction of published
scientific results in a traceable manner.  These are as listed

.. hlist::
    :columns: 4

    * `Arts 1.0.x <https://github.com/atmtools/arts/tree/v1.0.x>`_
    * `Arts 2.0.x <https://github.com/atmtools/arts/tree/v2.0.x>`_
    * `Arts 2.2.x <https://github.com/atmtools/arts/tree/v2.2.x>`_
    * `Arts 2.4.x (LTS) <https://github.com/atmtools/arts/tree/v2.4.x>`_
    * `Arts 2.4.0 <https://github.com/atmtools/arts/tree/v2.4.0>`_
    * `Arts 2.5.6 <https://github.com/atmtools/arts/tree/v2.5.6>`_
    * `Arts 2.5.8 <https://github.com/atmtools/arts/tree/v2.5.8>`_
    * `Arts 2.5.10 <https://github.com/atmtools/arts/tree/v2.5.10>`_

The API of Arts is constantly being updated to suit our and our users future need.
Do not expect code compatibility between any two patches.  Furthermore, only one LTS
version of Arts is ever supported in case bugs are found.
