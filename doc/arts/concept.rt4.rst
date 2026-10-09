.. _Sec RT4:

RT4 polarized doubling-adding solver
####################################

RT4 solves the polarized plane-parallel radiative-transfer equation with
thermal sources in an azimuthally symmetric medium, for the Stokes components
:math:`[I]` or :math:`[I,Q]`.  Unlike RT3 (:doc:`concept.rt3`), it takes the
optics of each layer as matrices on the streams, so they may describe
oriented particles with a polarized extinction.  It is RT3's doubling-adding
method in the azimuthal mode 0 :cite:p:`Evans1995`.  ARTS keeps RT4 as an
independent reference for VDISORT (:doc:`concept.disort`).

This page states the equations in the conventions of :doc:`concept.disort`
and :doc:`concept.rt3`, and the parts RT4 shares with RT3 are given there.
The implementation, interface and validation are documented in
:doc:`dev.rt4`.

Geometry and sign conventions
*****************************

The geometry, the streams and their weights, and the hemisphere notation are
those of RT3 (:ref:`Sec RT3 geometry`):
:math:`\mu>0` upwards, :math:`N` quadrature nodes :math:`\mu_j` with weights
:math:`w_j` integrating over :math:`[0,1]`, the same in both hemispheres,
and zero-weight extra directions after them.  Superscripts :math:`+` and
:math:`-` mark the upward and downward hemispheres; RT4's own "+" is
downward.  The layers are ordered top-down, and the vertical coordinate is
the height :math:`h`.

RT4 transports the leading :math:`n_s\le2` components of the ARTS Stokes
vector, :math:`[I]` or :math:`[I,Q]`, with :math:`Q=I_v-I_h` in the
meridional basis.

Polarized radiative transfer
****************************

Along a ray of direction cosine :math:`\mu`, :math:`\mathrm dh=\mu\,\mathrm ds`,
and the equation solved in a layer is

.. math::

   \mu\frac{\partial\boldsymbol I(h,\mu)}{\partial h}
   = -\boldsymbol K(\mu)\boldsymbol I
   + \boldsymbol a(\mu)B(T)
   + 2\pi\int_{-1}^{1}\overline{\boldsymbol Z}(\mu,\mu')\,
       \boldsymbol I(h,\mu')\,\mathrm d\mu',

with, per unit length,

* the extinction matrix :math:`\boldsymbol K(\mu)` of the particles plus the
  scalar gas extinction :math:`k_g\boldsymbol 1`,
* the absorption vector :math:`\boldsymbol a(\mu)` of the particles plus the
  gas absorption :math:`k_g\boldsymbol e_I`, and
* the azimuthal mean of the phase matrix, per unit length and steradian,

  .. math::

     \overline{\boldsymbol Z}(\mu,\mu')
     = \frac{1}{2\pi}\int_0^{2\pi}\boldsymbol Z(\mu,\mu',u)\,\mathrm du,

  which includes the number density of the particles.

:math:`B` is the Planck function of the temperature, linear in height within
a layer.  The field is azimuthally symmetric: there is no beam, and RT4
transports the azimuthal mean only, RT3's mode 0.

For a medium that is symmetric under reflection in a vertical plane (as for
randomly oriented particles, or particles oriented with random azimuths and a
plane of symmetry), :math:`\boldsymbol K=\boldsymbol D\boldsymbol K\boldsymbol D`
with RT3's :math:`\boldsymbol D=\operatorname{diag}(1,1,-1,-1)`, so the
blocks of :math:`\boldsymbol K` that couple :math:`(I,Q)` to :math:`(U,V)`
vanish, and those blocks of :math:`\boldsymbol Z` are odd in the azimuth, so
their means in :math:`\overline{\boldsymbol Z}` vanish.  :math:`U` and
:math:`V` then decouple from :math:`I` and :math:`Q` and have no source, so
they are zero, and :math:`n_s=2` is exact.

Stream equations
================

On the quadrature, with the signed streams :math:`\mu_i` of both hemispheres
and the repeated weights :math:`\overline w_j`,

.. math::

   \mu_i\frac{\mathrm d\boldsymbol I_i}{\mathrm dh}
   = -\boldsymbol K(\mu_i)\boldsymbol I_i
   + \boldsymbol a(\mu_i)B
   + 2\pi\sum_j\overline w_j\,\overline{\boldsymbol Z}(\mu_i,\mu_j)\boldsymbol I_j.

Energy is conserved on the streams only if every incident quadrature stream
scatters the difference of extinction and absorption into the quadrature
streams of both hemispheres:

.. math::

   K_{11}(\pm\mu_j) = a_1(\pm\mu_j)
   + 2\pi\sum_{i=1}^{N}w_i
     \left[\overline{\boldsymbol Z}(\mu_i,\pm\mu_j)
          +\overline{\boldsymbol Z}(-\mu_i,\pm\mu_j)\right]_{II}.

``rt4::solve`` rejects optics that miss it by more than
``normalisation_tolerance`` times :math:`K_{11}`.

Doubling and adding
*******************

A slab is represented by its reflection, transmission and source operators,
and layers are built by doubling and combined by adding exactly as in RT3
(:ref:`Sec RT3 doubling-adding`), with the
thermal source only.  Only the initial sublayer differs.

Initial layer
=============

A scattering layer of thickness :math:`\Delta h` is halved :math:`n_d` times
as in RT3, with the vertical optical thickness
:math:`\Delta\tau=(K_{11}(-\mu_1)+k_g)\Delta h` taken from the first
downward stream.  To first order in the sublayer thickness :math:`\delta h`,

.. math::

   \boldsymbol T^{\pm\pm}_{ij}
   &= \delta_{ij}\boldsymbol 1
    - \frac{\delta h}{\mu_i}\left(\delta_{ij}\boldsymbol K(\pm\mu_i)
    - 2\pi w_j\,\overline{\boldsymbol Z}(\pm\mu_i,\pm\mu_j)\right),\\
   \boldsymbol R^{\pm\mp}_{ij}
   &= \frac{\delta h}{\mu_i}\,2\pi w_j\,\overline{\boldsymbol Z}(\pm\mu_i,\mp\mu_j),\\
   \boldsymbol S^{\pm}_i
   &= \frac{\delta h}{\mu_i}\,\boldsymbol a(\pm\mu_i)\,B_{\mathrm t},

with :math:`\mu_i>0` the cosine of the outgoing stream, :math:`\boldsymbol K`
and :math:`\boldsymbol a` including the gas, and :math:`B_{\mathrm t}` the
Planck function at the top of the layer.  The doubling raises the Planck
function across the sublayers as RT3 does, a staircase on the linear profile.

Mirror symmetry between the hemispheres
=======================================

RT4's doubling computes only one direction of each slab: it assumes
:math:`\boldsymbol T^{++}=\boldsymbol T^{--}` and
:math:`\boldsymbol R^{+-}=\boldsymbol R^{-+}`.  That holds when the medium is
mirror symmetric between the hemispheres,

.. math::

   \boldsymbol K(-\mu) = \boldsymbol K(\mu),\qquad
   \boldsymbol a(-\mu) = \boldsymbol a(\mu),\qquad
   \overline{\boldsymbol Z}(-\mu_o,-\mu_i) = \overline{\boldsymbol Z}(\mu_o,\mu_i),

the last for both the transmitted (:math:`\mu_o\mu_i>0`) and the reflected
(:math:`\mu_o\mu_i<0`) parts.  For :math:`n_s\le2` this is RT3's
:math:`\boldsymbol D`-symmetry with :math:`\boldsymbol D` the identity.
``rt4::solve`` rejects optics that break it.

Absorbing layers
================

A layer without particles is integrated exactly, with the formulas of RT3
(:ref:`Sec RT3 absorbing layers`) for the optical
thickness :math:`k_g\Delta h`.

Sources and boundary conditions
*******************************

* The only internal source is the thermal emission
  :math:`\boldsymbol a B`: unpolarized for the gas, polarized as
  :math:`\boldsymbol a` for the particles.
* The radiance incident at the top is
  :math:`\boldsymbol I^-_{\mathrm{top}}=B(T_{\mathrm{sky}})\boldsymbol e_I`.

Surface reflection
******************

The surface is a slab of its own as in RT3: at the surface
:math:`\boldsymbol I^+=\boldsymbol R_{\mathrm s}\boldsymbol I^-+\boldsymbol I_{\mathrm g}`,
with

* **Lambertian** (albedo :math:`A`):
  :math:`(\boldsymbol R_{\mathrm s})_{ij}=2A\,\mu_jw_j\,\boldsymbol e_I\boldsymbol e_I^{\mathsf T}`,
  :math:`\boldsymbol I_{\mathrm g}=(1-A)B(T_{\mathrm s})\boldsymbol e_I`;
* **Fresnel**: RT3's :math:`\boldsymbol R_F(\mu_i)` restricted to
  :math:`[I,Q]`, specular, with
  :math:`\boldsymbol I_{\mathrm g}=[(1-R_1)B,\,-R_2B]^{\mathsf T}`;
* **Specular** with a fixed reflection matrix :math:`\boldsymbol R`:
  :math:`(\boldsymbol R_{\mathrm s})_{ij}=\delta_{ij}\boldsymbol R`,
  :math:`\boldsymbol I_{\mathrm g}=(\boldsymbol 1-\boldsymbol R)B(T_{\mathrm s})\boldsymbol e_I`;
* **Discrete**: :math:`\boldsymbol R_{\mathrm s}` and
  :math:`\boldsymbol I_{\mathrm g}` given directly on the streams, quadrature
  factors included.

Radiances at the levels
***********************

The radiances at every level follow from the slabs above and below it as in
RT3 (:ref:`Sec RT3 levels`).
They are those of the azimuthal mean; RT4 returns no fluxes.

Mathematical limitations
************************

* **First-order initial layer.**  As in RT3, the error of a scattering layer
  is first order in the sublayer thickness, and the doublings are chosen
  from the vertical optical thickness at the first stream, so the slant
  thickness at small :math:`\mu_i` is larger.
* **Azimuthal symmetry.**  There is no beam, and only :math:`[I]` or
  :math:`[I,Q]` are transported, exact for media symmetric under reflection
  in a vertical plane.
* **Mirror symmetry between the hemispheres** is required by the doubling.
* **Normalization on the streams** must hold for the given optics; phase
  matrices with forward peaks between the streams need more streams or
  renormalization before RT4.
* **Plane-parallel**, as RT3.

See :doc:`dev.rt4` for the implementation and its validation.
