.. _Sec Scattering species:

Scattering species
##################

A scattering species is a population of scatterers whose bulk scattering
properties, the extinction, the absorption and the phase matrix, follow from
the local state of the medium it is in.  The workspace holds two arrays of
them, :class:`~pyarts3.arts.ArrayOfScatteringSpecies`:

- :attr:`~pyarts3.workspace.Workspace.atm_scat_species` for the atmosphere,
- :attr:`~pyarts3.workspace.Workspace.subsurf_scat_species` for the subsurface.

The bulk scattering properties of an array are the sums over its species.

Points and spectral points
**************************

A species is evaluated at a point of its medium together with the spectral
point at the same position:

- in the atmosphere, an :class:`~pyarts3.arts.AtmPoint` and an
  :class:`~pyarts3.arts.AtmSpectralPoint`;
- in the subsurface, a :class:`~pyarts3.arts.SubsurfacePoint` and a
  :class:`~pyarts3.arts.SubsurfaceSpectralPoint`.

The point holds the local state of the medium, see :ref:`Sec Atmospheric field`
and :ref:`Sec Subsurface field`.  The spectral point holds the properties of
the spectral field, :attr:`~pyarts3.workspace.Workspace.atm_spectral_field` or
:attr:`~pyarts3.workspace.Workspace.subsurf_spectral_field`, on the frequencies
of the calculation.  A method that evaluates the species takes both at every
position it uses.  Each species reads only what it needs, by the names it holds
or by the keys of its model, and a species that needs what its medium does not
have raises an error.

Types of scattering species
***************************

- :class:`~pyarts3.arts.HenyeyGreensteinScatterer` holds the names of three
  spectral properties: the extinction [m\ :sup:`-1`], the single scattering
  albedo and the asymmetry parameter :math:`g`.  Its phase function is that
  of Henyey and Greenstein, :math:`p`, and its scattering matrix has
  :math:`F_{11}=F_{22}=p`,
  :math:`F_{33}=F_{44}=p\,(3\cos\Theta-\cos^3\Theta)/2` and
  :math:`F_{12}=F_{34}=0`: it produces no polarization and is regular at
  forward and backward scattering.  It works in both media.
- :class:`~pyarts3.arts.GasScatterer` is molecular scattering, for example
  :func:`~pyarts3.arts.GasScatterer.air_simple_rayleigh`.  It takes the number
  density of the gas from the pressure and temperature of the atmospheric
  point, so it works in the atmosphere only.
- :class:`~pyarts3.arts.ScatteringHabit` combines the single scattering data
  of a :class:`~pyarts3.arts.ParticleHabit`, for example from the T-matrix
  method (:doc:`user.tmatrix`), with a particle size distribution.  The
  parameters of the distribution are properties of the atmospheric field,
  identified by :class:`~pyarts3.arts.ScatteringSpeciesProperty`, so it works
  in the atmosphere only.
- :class:`~pyarts3.arts.ScatteringGeneralSpectralTRO` is a Python function of
  the point, the spectral point, the frequency grid and the degree, which
  returns the bulk scattering properties as a Legendre series,
  :class:`~pyarts3.arts.ScatteringTroSpectralVector`.  The point is an
  :class:`~pyarts3.arts.AtmPoint` or a :class:`~pyarts3.arts.SubsurfacePoint`,
  depending on where the species is used, and the function decides which it
  supports.  Its ``degree`` is that of the series where the calling method
  does not choose it.

Here an atmospheric ice species is a Henyey-Greenstein scatterer, and its
properties are evaluated at two frequencies:

.. code-block:: python

  import pyarts3 as pyarts

  ice = pyarts.arts.HenyeyGreensteinScatterer("ice extinction", "ice ssa", "ice g")
  atm_scat_species = pyarts.arts.ArrayOfScatteringSpecies()
  atm_scat_species.add(ice)

  atm_spectral_point = pyarts.arts.AtmSpectralPoint()
  atm_spectral_point["ice extinction"] = [1e-4, 2e-4]
  atm_spectral_point["ice ssa"] = [0.9, 0.9]
  atm_spectral_point["ice g"] = [0.5, 0.7]

  bulk = atm_scat_species.tro_spectral(pyarts.arts.AtmPoint(), atm_spectral_point, [89e9, 183e9], 8)
  bulk.extinction_matrix  # the extinction of each frequency

In a workspace, the same names are properties of
:attr:`~pyarts3.workspace.Workspace.atm_spectral_field`, see
:ref:`Sec Atmospheric field`, and the methods evaluate the spectral points
themselves.  The same species as a Python function:

.. code-block:: python

  def ice_series(point, spectral_point, f_grid, degree):
      return ice.tro_spectral(point, spectral_point, f_grid, degree)

  atm_scat_species.add(pyarts.arts.ScatteringGeneralSpectralTRO(ice_series, 8))

Where the scattering species are used
*************************************

:attr:`~pyarts3.workspace.Workspace.atm_scat_species` is used by the
scattering of the propagation matrices along a path
(``spectral_propmat_scatAddSpectralScatteringSpeciesTRO``), and through them
by the DISORT and VDISORT settings with the "ScatteringSpecies" option, by
first-order solar scattering, by Monte Carlo and by the radar methods.
:attr:`~pyarts3.workspace.Workspace.subsurf_scat_species` is the volume
scattering of the subsurface for VDISORT, see :doc:`user.disort`.

Derivatives
***********

A species gives the derivatives of its bulk scattering properties by the
quantities it reads: a Henyey-Greenstein scatterer by its three spectral
properties, and a habit by the temperature and the properties of its size
distribution.  The derivatives by everything it does not read are zero.
Retrieval targets on the spectral fields are described in :doc:`user.oem`.
