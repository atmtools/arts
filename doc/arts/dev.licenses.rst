Licences of bundled code
========================

ARTS compiles in code from other projects, mostly from ``3rdparty/``.  Every
distribution of ARTS, from a source tree to a precompiled package, must carry
the licence notices of the bundled code it contains.  ARTS itself is
LGPL-3.0-or-later OR GPL-3.0-or-later (``LICENSE.txt``).

Adding bundled code
-------------------

#. Keep the component's licence and notice files unchanged next to its
   code.  If its copyright and permission notice exist only in source
   headers, add a verbatim copy as ``LICENSE`` or ``NOTICE`` in its
   directory.
#. Register it in the CMake branch that builds or links it, so that exactly
   the components of the current configuration are registered:

   .. code-block:: cmake

      arts_add_bundled_license(polradtran SPDX MIT FILES LICENSE)

   ``SPDX`` is the component's SPDX licence expression, or
   ``LicenseRef-<id>`` for a licence without an SPDX identifier.  ``FILES``
   are relative to the calling directory.  The function is defined in
   ``cmake/modules/ArtsBundledLicenses.cmake``.

   The registered files are configure dependencies of the build.  If one is
   missing, configuring fails with ``The licence file ... of the bundled
   component ... does not exist``; if one is deleted after configuring, the
   next build reruns CMake and fails the same way.  Editing one reconfigures
   and recompiles its compiled-in copy.  Registering the same name twice also
   fails.
#. Add the component to the table below.

What a build carries
--------------------

The registered components of a configuration are gathered into the build:

* **Compiled into ARTS.**  ``src/core/licenses`` generates a source with the
  complete licence files of every registered component.  Every binary
  carries them, and they are available from C++ (``bundled_licenses.h``)
  and Python:

  .. code-block:: python

     import pyarts3 as pyarts

     print(pyarts.arts.globals.license_expression())
     for component in pyarts.arts.globals.bundled_components():
         print(component.name, component.spdx, list(component.files))

* **In the Python package.**  The package contains every registered file as
  ``licenses/<name>/<file>`` with the index ``licenses/THIRD_PARTY.txt``,
  next to ``LICENSE.txt``.  Its ``License-Expression`` metadata is the
  licence of ARTS AND those of the components, for example
  ``(LGPL-3.0-or-later OR GPL-3.0-or-later) AND Apache-2.0 WITH
  LLVM-exception AND ... AND MIT``.  Configuring prints it as
  ``Python package licence``.

``tests/core/licenses/bundled-licenses.py`` checks that the compiled-in
components are the ones this build contains, and that they match the
package's files.

Which licences go into which builds
-----------------------------------

* Permissive licences (MIT, BSD, Apache-2.0) and LGPL-3.0-or-later can be in
  every build, including ``ENABLE_ARTS_LGPL`` builds.  Such components may be
  built by default.
* GPL and other strong copyleft licences can only be in builds without
  ``ENABLE_ARTS_LGPL``, which turns those components off.
* A component whose licence is not an open-source licence is opt-in only.

Current components:

.. list-table::
   :header-rows: 1

   * - Component
     - Licence (SPDX)
     - Files
     - Compiled in
   * - invlib
     - MIT
     - ``LICENSE.txt``
     - always
   * - mdspan
     - Apache-2.0 WITH LLVM-exception
     - ``LICENSE``
     - always
   * - Faddeeva
     - MIT
     - ``LICENSE``
     - always
   * - wigxjpf
     - LGPL-3.0-or-later
     - ``NOTICE``, ``COPYING.LESSER``, ``COPYING``
     - always
   * - fastwigxj
     - LGPL-3.0-or-later
     - ``NOTICE``, ``COPYING.LESSER``, ``COPYING``
     - with ``FASTWIGNER``
   * - polradtran (RT3, RT4)
     - MIT
     - ``LICENSE``
     - when a Fortran compiler is found, unless ``ENABLE_RT3`` and
       ``ENABLE_RT4`` are off
   * - cdisort
     - GPL-3.0-or-later
     - ``NOTICE``, ``COPYING``
     - with ``ENABLE_CDISORT``, not in LGPL builds
   * - shtns
     - CeCILL-2.1
     - ``COPYRIGHT``, ``LICENSE``
     - when SHTns is used, not in LGPL builds
   * - tmatrix
     - LicenseRef-Mishchenko-T-matrix (free use in not-for-profit research)
     - ``License.txt``
     - opt-in with ``ENABLE_TMATRIX``, not in LGPL builds

The notice in cdisort's ``cdisort.h`` grants GPL version 3 or later, while
its ``COPYING`` holds the GPL version 2 text; ``NOTICE`` copies the header,
and the GPL version 3 text is in ARTS's ``LICENSE.txt``.

External libraries
------------------

Libraries that ARTS links dynamically are not bundled code and are not
registered: the BLAS and LAPACK implementation, FFTW (GPL-2.0-or-later, when
SHTns is used), the Fortran runtime (libgfortran and libquadmath,
GPL-3.0-or-later WITH GCC-exception-3.1, when RT3, RT4 or T-matrix is built),
and the C++ and OpenMP runtimes.  A conda package declares them as
dependencies, and their packages carry their own licences.  A tool that
copies shared libraries into a package, such as ``auditwheel``,
``delocate`` or ``delvewheel``, does not add their licences, so a package
built that way must ship those licences as well.
