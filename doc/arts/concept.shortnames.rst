Short-hand names for ARTS variables
###################################

All variables in ARTS use short-hand names to identify them and
simplify writing. The keys to read these short-hand names are listed
below.

All complete short-hand names are separated by an underscore (``_``).

.. rubric:: Key components

A combination of these components form the main part of the variable names.

.. list-table::
   :header-rows: 1

   * - Short name
     - Description
   * - ``atm``
     - Used for atmospheric variables.
   * - ``surf``
     - Used for surface variables.
   * - ``subsurf``
     - Used for subsurface variables.

   * - ``abs``
     - Used for absorption-related variables.
   * - ``rad``
     - Used for radiance-related variables.
   * - ``propmat``
     - Used for propagation matrix variables.
   * - ``tranmat``
     - Used for transmission matrix variables.
   * - ``phamat``
     - Used for phase matrix variables.
   * - ``srcvec``
     - Used for source vector variables.
   * - ``absvec``
     - Used for absorption vector variables.
   * - ``scat``
     - Used for scattering-related variables.
   * - ``refl``
     - Used for reflection-related variables.

   * - ``jac``
     - Used for Jacobian variables.
   * - ``covmat``
     - Used for covariance-matrix variables.
   * - ``leg``
     - Used for Legendre variables.

   * - ``nlte``
     - Used for non-LTE-related variables.
   * - ``cia``
     - Used for collision-induced absorption-related variables.
   * - ``ecs``
     - Used for error-corrected sudden-related variables.
   * - ``xfit``
     - Used for cross-section fitting-related variables.
   * - ``meas``
     - Used for measurement-related variables.

   * - ``alt``
     - Used for altitude-related variables.
   * - ``lat``
     - Used for latitude-related variables.
   * - ``lon``
     - Used for longitude-related variables.
   * - ``freq``
     - Used for frequency-related variables.
   * - ``za``
     - Used for zenith angle-related variables.
   * - ``az``
     - Used for azimuth angle-related variables.
   * - ``ray``
     - Used for ray tracing-related variables.
   * - ``bkg``
     - Used for background-related variables.
   * - ``obs``
     - Used for observation-related variables.
   * - ``pos``
     - Used for position-related variables.
   * - ``los``
     - Used for line-of-sight-related variables.

.. rubric:: Structural components

These are generally the last or first part of the variable names
if they are present.

.. list-table::
   :header-rows: 1

   * - Short name
     - Description
     - Position
   * - ``grid``
     - Used for grid variables.  Grid variables define the grids on which
       other variables are defined.  Grids may be spatial grids (e.g.,
       altitude grid) or other types of grids (e.g., frequency grid).
     - At the end.
   * - ``field``
     - Used for field variables.  Field variables are defined on
       one or more grids, i.e., they are multi-dimensional arrays.
     - At the end.
   * - ``point``
     - Used for point variables.  Point variables are single values.
       You can imagine them as extractions from the field variables at
       specific grid points.  Note that grid variables are not considered
       point variables, e.g., an altitude grid does not consist of altitude points.
     - At the end.
   * - ``profile``
     - Used for profile variables.  Profiles are one-dimensional arrays.
       They are composed of point values orderly extracted on a specific
       grid.  Otherwise they are similar to field variables.
     - At the end.
   * - ``path``
     - Used for variables defined along a path.  Path variables are like
       profile variables, but they are not defined along a single grid axis.
       Instead, they are defined along a trajectory in the grid space.
     - At the end.
   * - ``spec``
     - Used for variables defined along a spectral axis.  That is, along
       the frequency grid.
     - At the beginning.
   * - ``single``
     - Used for variables defined along a spectral axis, however only a single
       frequency point is being considered.
     - At the beginning.
