Licences of bundled code
========================

ARTS compiles in code from other projects, mostly from ``3rdparty/``.  Every
distribution of ARTS, from a source tree to a precompiled package, must carry
the licence notices of the bundled code it contains.  ARTS itself is
LGPL-3.0-or-later OR GPL-3.0-or-later (``LICENSE.txt``).

The mechanism is in ``cmake/modules/ArtsBundledLicenses.cmake``, whose header
documents its functions and their arguments.

Adding bundled code
-------------------

#. Keep the component's licence and notice files unchanged next to its
   code.  If its copyright and permission notice exist only in source
   headers, add a verbatim copy as a ``LICENSE`` or ``NOTICE`` file in its
   directory.
#. Register the component, with its SPDX licence expression and its licence
   files, in the branch of the CMake logic that builds or links it.  Exactly
   the components of the current configuration are then registered.  A
   licence without an SPDX identifier gets a ``LicenseRef-`` identifier.
   The existing registrations in the ``CMakeLists.txt`` files of
   ``3rdparty/`` show the pattern.
#. The registered files are configure dependencies:

   * A missing file fails the configuration, and so does a name registered
     twice.
   * A file deleted after configuring makes the next build rerun CMake and
     fail.
   * Editing a file reconfigures, and recompiles its compiled-in copy.

What a build carries
--------------------

The registered components of a configuration are gathered into the build:

* **Compiled into ARTS.**  ``src/core/licenses`` generates a source with the
  complete licence files of every registered component, so every binary
  carries them.  They are available from C++ (``bundled_licenses.h``) and
  from Python, in the ``globals`` submodule of ``pyarts3.arts``, together
  with the licence expression of the build.
* **In the Python package.**  The package contains every registered file
  under ``licenses/``, one directory per component, with an index
  ``licenses/THIRD_PARTY.txt``, next to ``LICENSE.txt``.  Its
  ``License-Expression`` metadata is the licence of ARTS AND those of the
  components, and configuring prints it.

The Python test in ``tests/core/licenses`` checks that the compiled-in
components are the ones this build contains, and that they match the
package's files.  A build's own list of components, with their licences, is
its ``THIRD_PARTY.txt`` or the Python listing.  The registrations are the
source.

Which licences go into which builds
-----------------------------------

* Permissive licences (MIT, BSD, Apache-2.0) and LGPL-3.0-or-later can be in
  every build, including ``ENABLE_ARTS_LGPL`` builds.  Such components may be
  built by default.
* GPL and other strong copyleft licences can only be in builds without
  ``ENABLE_ARTS_LGPL``, which turns those components off.
* A component whose licence is not an open-source licence is opt-in only.

A component whose notice and shipped licence text disagree keeps both.  The
notice is copied verbatim into a ``NOTICE`` file, and the registered SPDX
expression follows the notice.  A comment at the registration explains the
discrepancy.

External libraries
------------------

Libraries that ARTS links dynamically are not bundled code, and they are not
registered: the BLAS and LAPACK implementation, FFTW, the Fortran runtime
when Fortran code is built, and the C++ and OpenMP runtimes.  A conda package
declares them as dependencies, and their packages carry their own licences.
A tool that copies shared libraries into a package, such as ``auditwheel``,
``delocate`` or ``delvewheel``, does not add their licences, so a package
built that way must ship those licences as well.
