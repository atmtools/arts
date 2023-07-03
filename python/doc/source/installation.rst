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

.. list-table:: Versions with python support similar to today
    :header-rows: 0
    :widths: 25 75

    * - `Arts 2.5.10 <https://github.com/atmtools/arts/tree/v2.5.10>`_
      - Latest stable version
    * - `Arts 2.5.8 <https://github.com/atmtools/arts/tree/v2.5.8>`_
      -
    * - `Arts 2.5.6 <https://github.com/atmtools/arts/tree/v2.5.6>`_
      -

.. list-table:: Before current python support
    :header-rows: 0
    :widths: 25 75

    * - `Arts 2.4.x <https://github.com/atmtools/arts/tree/v2.4.x>`_
      -
    * - `Arts 2.4.0 <https://github.com/atmtools/arts/tree/v2.4.0>`_
      -
    * - `Arts 2.2.x <https://github.com/atmtools/arts/tree/v2.2.x>`_
      - See `planetary toobox paper <https://doi.org/10.5194/gmd-11-1537-2018>`_
    * - `Arts 2.0.x <https://github.com/atmtools/arts/tree/v2.0.x>`_
      - See `version 2 paper <https://doi.org/10.1016/j.jqsrt.2011.03.001>`_
    * - `Arts 1.0.x <https://github.com/atmtools/arts/tree/v1.0.x>`_
      - See `version 1 paper <https://doi.org/10.1016/j.jqsrt.2004.05.051>`_

The API of Arts is constantly being updated to suit our and our users future need.
Do not expect code compatibility between any two patches.  We strongly recommend to
use the latest version of Arts on conda unless you want to help with the development.
