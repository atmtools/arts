.. _Sec Surface scattering concept:

Surface scattering framework
############################

This page develops the mathematics and physics of the ARTS surface scattering
framework, the BRDF-based lower boundary condition of the radiative transfer
equation.  The user interface is described in :doc:`user.surface_scattering`,
the implementation in :doc:`dev.surface_scattering`.

The boundary condition
**********************

When a ray path terminates at the surface, the upwelling spectral radiance is
the sum of the incoming radiance scattered by the surface and the surface or
sub-surface emission attenuated by the surface's own transparency:

.. math::

  \mathbf{I}_\text{out}(\nu, \Omega_\text{out})
    = \sum_{\Omega_\text{in}}
        \mathbf{BRDF}(\nu, \Omega_\text{in}, \Omega_\text{out})
        \mathbf{I}_\text{in}(\nu, \Omega_\text{in})
    + \boldsymbol{\varepsilon}(\nu, \Omega_\text{out})
      \mathbf{I}_\text{surface}(\nu, \Omega_\text{out})

All quantities are fully polarimetric: :math:`\mathbf{I}` are Stokes vectors,
while the BRDF and the emissivity :math:`\boldsymbol{\varepsilon}` are 4×4
Mueller tensors.  All angles are measured with respect to the *actual* local
surface normal supplied by the surface field, so tilted surfaces are handled
exactly, not approximated as flat.

The sum over :math:`\Omega_\text{in}` is a discretised integral over the
incoming hemisphere; how it is sampled is a design choice of the framework and
is discussed in `Sampling of the incoming radiation`_ below.

Model families and Kirchhoff consistency
****************************************

A surface scattering model provides a BRDF and, complementarily, an emissivity.
Every built-in model ships the *complementary* emissivity built in, so
Kirchhoff's law :math:`\mathbf{BRDF} + \boldsymbol{\varepsilon} = \mathbf{I}_4`
is enforced at the model level, not by the radiative transfer machinery.  The
three families are:

Lambertian (diffuse)
  Isotropic scattering with spectrally varying reflectivity :math:`r(\nu)`:

  .. math::

    \mathbf{BRDF}_\text{diffuse} = \mathrm{diag}(r, 0, 0, 0),
    \qquad
    \boldsymbol{\varepsilon}_\text{diffuse} = \mathrm{diag}(1-r, 0, 0, 0)

  Only the intensity channel is scattered or emitted; the polarisation
  channels are untouched.

Fresnel (specular)
  Specular reflection from a medium of real refractive index :math:`n_2(\nu)`
  under an incident medium of index :math:`n_1`, computed from the Fresnel
  amplitude reflectances :math:`R_v, R_h` at the incidence angle:

  .. math::

    \mathbf{BRDF}_\text{specular} = \mathbf{M}(R_v, R_h),
    \qquad
    \boldsymbol{\varepsilon}_\text{specular} = \mathbf{I}_4 - \mathbf{M}(R_v, R_h)

  This is the polarimetric Fresnel reflectance Mueller matrix (the same tensor
  as in :meth:`~pyarts3.workspace.Workspace.spectral_surf_reflFlatRealFresnel`).

Flat scalar (specular)
  Angle-independent scalar reflectance :math:`R(\nu) \in [0, 1]` applied
  identically to all Mueller channels:

  .. math::

    \mathbf{BRDF}_\text{specular} = \mathbf{M}\{R\},
    \qquad
    \boldsymbol{\varepsilon}_\text{specular} = \mathbf{M}\{1 - R\}

Each model fills only the tensor pair of its own family — a diffuse model
leaves the specular pair empty and vice versa — so a Lambertian and a specular
model can be blended at the same point without interference: their
contributions simply add in different tensor channels.

Units convention
****************

The stored BRDFs are *dimensionless scattering kernels*, not physical BRDFs
in :math:`[\mathrm{sr}^{-1}]`: the physical Lambertian BRDF would be
:math:`r/\pi`, but ARTS stores just :math:`r`.  Consistently, the angular
quadrature weights are *raw multipliers* — the scattering methods insert no
geometric factors themselves.  For physical flux conservation the weights must
therefore carry the projected solid angle factor :math:`\cos(\theta_\text{in})
\,\mathrm{d}\Omega`; with plain solid-angle weights :math:`\sin(\theta)\,
\mathrm{d}\theta\,\mathrm{d}\phi` a uniform incoming radiance field would
reflect with effective reflectance :math:`2\pi r` instead of :math:`r`.

The two exceptions are built in: the sun-beam diffuse method folds the
projected-area factor :math:`\cos(\theta_\text{inc})` in itself (a single sun
beam scatters as :math:`r \cos(\theta_\text{inc}) \mathbf{I}_\text{beam}`),
and the sun-beam specular method needs no weight at all because specular
reflection conserves radiance along the mirror direction.

Sampling of the incoming radiation
**********************************

The four scattering methods of the framework differ in *which* tensor pair they
contract and *how* they sample the incoming direction integral:

.. list-table::
   :header-rows: 1

   * - 
     - Diffuse tensors
     - Specular tensors
   * - **Hemispheric sampling**
     - Quadrature over a user grid :math:`\theta \times \phi`
     - Single mirror direction of the outgoing ray
   * - **Sun beams**
     - One delta beam per sun, with :math:`\cos\theta_\text{inc}` folded in
     - Contributes only if the mirror direction hits the solar disc

Two geometric facts govern all four:

*Horizon gating.*  Incoming directions are checked against the horizon defined
by the *actual* surface normal.  Directions below that horizon contribute
nothing; directions that a flat horizon would hide but the tilt exposes *are*
included.  For the quadrature this means the zenith grid must extend beyond
90° to capture directions made visible by a surface tilt.

*The solar disc.*  The sun is an extended source of angular radius
:math:`\alpha = \arcsin(R_\text{sun}/d_\text{sun}) \approx 0.27°`.  A practical
quadrature cannot resolve it, so suns enter exactly only through the sun-beam
methods, which test whether the relevant direction falls inside the disc —
first geometrically, then on the refraction-aware beam line-of-sight obtained
from the observer path agenda.  The same disc test gates the local horizon for
beams: a sun below the tilted horizon contributes nothing.

Additive structure and double counting
**************************************

The scattering methods are *additive*: each adds its contribution to
``spectral_rad`` and ``spectral_rad_jac``, after an initialisation method sizes
and zeroes them.  This makes the boundary condition composable — any subset of
the four methods can be chained in any order behind one initialisation.

The price of composition is two bookkeeping obligations, because the emission
term and the sun are reachable by more than one method:

*Emission.*  Each method adds its tensor family's emission term.  Pairing a
hemispheric method with its sun-beam counterpart must therefore disable the
emission of the sun-beam method, or the emission is counted twice.

*The sun.*  With suns present, the quadrature or mirror direction may also
contain the solar disc.  When a sun-beam method is chained, the hemispheric
methods must *exclude* directions containing a sun, or the sun is counted
twice.  Without a sun-beam method in the chain, the exclusion must stay off —
the hemispheric route is then the only route to the sun.

With both obligations honoured, the chain satisfies the exact identity

.. math::

  \mathbf{I}_\text{full}
    = \mathbf{I}_\text{hemispheric}
    + \mathbf{I}_\text{sun-beams}
    - \boldsymbol{\varepsilon}\,\mathbf{I}_\text{surface}

to machine precision, whenever no traced direction already contains a sun.

Jacobian propagation
********************

The Jacobian :math:`\partial \mathbf{I}_\text{out} / \partial \mathbf{x}`
propagates through exactly two routes, with the surface tensors acting as
constant coefficients:

.. math::

  \frac{\partial \mathbf{I}_\text{out}}{\partial \mathbf{x}}
    = \sum_{\Omega_\text{in}} \mathbf{BRDF}\,
        \frac{\partial \mathbf{I}_\text{in}}{\partial \mathbf{x}}
    + \boldsymbol{\varepsilon}\,
      \frac{\partial \mathbf{I}_\text{surface}}{\partial \mathbf{x}}

The first term comes from the agenda computing the incoming radiance (the
atmospheric state, and anything else along the incoming paths); the second
from the agenda computing the surface emission (surface temperature,
sub-surface scattering properties).

.. note::

   There is **no Jacobian with respect to the scattering model itself**:
   reflectivities, refractive indices, and the surface mask weights are not
   differentiable parameters.  In retrievals they act as fixed prior
   information.
