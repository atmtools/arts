.. _Workspace Dimensions:

Workspace Dimensions
====================

A dimension is a size that several workspace variables share.  Naming one says
that the variables agree about it: every variable with ``NFREQ`` has as many
frequency points as every other local variable, so a propagation matrix computed on one
frequency grid cannot be combined with a radiance computed on another.

Each workspace variable that has a shape reports it as an *effective shape*,
which lists its dimensions from the outside in.  For an array whose elements
are all shaped alike, the effective shape covers both: *spectral_propmat_path*
is ``[NPATH, NFREQ]`` because it holds one propagation matrix vector per path point,
each of them over the same frequency grid.  So the effective shape is what an
index into the variable has to provide.

These are not free-form annotations.  Methods and agendas verify them when a
user supplies the data or composes the calls, and the failure says which
variables disagreed and what sizes they had.  Every method and agenda that
verifies sizes lists what it requires under its own *Constraints* rubric.

Nothing here says how a size is read.  That belongs to the group of the
variable, since a size is read differently from a :class:`~pyarts3.arts.Vector`
than from a :class:`~pyarts3.arts.JacobianTargets`, while the dimension is the
same in both.

.. include:: workspace.dimensions.auto.rst
