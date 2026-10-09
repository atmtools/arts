.. _Sec RT3:

RT3 polarized doubling-adding solver
####################################

RT3 solves the polarized plane-parallel radiative-transfer equation with a
solar beam and thermal sources for randomly oriented particles with a plane
of symmetry.  It separates the field into azimuthal Fourier modes, and for
every mode it builds the reflection, transmission and source operators of
each layer by doubling a thin initial layer, then combines the layers by
adding.  The method is that of :cite:t:`Evans1999`.  ARTS keeps RT3 as an
independent reference for VDISORT (:doc:`concept.disort`); its RT4 sibling,
restricted to azimuthally symmetric problems, is described in
:doc:`concept.rt4`.

This page states the equations in the conventions used for DISORT and VDISORT
in ARTS.  The implementation, interface and validation are documented in
:doc:`dev.rt3`.

.. _Sec RT3 geometry:

Geometry and sign conventions
*****************************

The atmosphere consists of homogeneous plane-parallel layers.  Optical depth
:math:`\tau` is zero at the top and increases downwards.  Direction cosines
are defined by

.. math::

   \mu = \cos\theta,

with :math:`\mu>0` directed upwards and :math:`\mu<0` directed downwards, as
in :doc:`concept.disort`.  The direct beam propagates downwards with direction
cosine :math:`-\mu_0`, :math:`\mu_0>0`, towards the azimuth :math:`\phi=0`.
The azimuth :math:`\phi` of a ray is that of its propagation direction,
counted counterclockwise seen from above.  ARTS's propagation directions
:math:`(\theta_z,\alpha)`, with the azimuth :math:`\alpha` clockwise seen from
above, are therefore

.. math::

   \mu = \cos\theta_z,\qquad \phi = -\alpha,

and VDISORT's radiance at the azimuth :math:`\phi_0+\phi` is RT3's at
:math:`\phi`.

RT3 transports the ARTS Stokes vector

.. math::

   \boldsymbol I = [I,Q,U,V]^{\mathsf T},\qquad Q = I_v - I_h,

in the meridional basis of each ray (the vertical polarization lies in the
plane of the ray and the zenith), the same in both hemispheres.  The leading
:math:`n_s` components are transported, :math:`n_s=1,\ldots,4`.

The streams are :math:`N` quadrature nodes :math:`\mu_1,\ldots,\mu_N` in
:math:`(0,1]` with weights :math:`w_j`, used in both hemispheres.  The weights
integrate over :math:`[0,1]`,

.. math::

   \sum_{j=1}^N w_j g(\mu_j) \approx \int_0^1 g(\mu)\,\mathrm d\mu,\qquad
   \sum_{j=1}^N w_j = 1.

The rules are the positive half of a :math:`2N`-point Gauss--Legendre rule
(Gauss), an :math:`N`-point Gauss--Legendre rule on :math:`[0,1]` (double
Gauss), or the positive half of a :math:`2N`-point Lobatto rule (Lobatto).
Extra directions with zero weight may follow the nodes; see
`Arbitrary-angle radiances`_.

Below, a superscript :math:`+` marks the upward hemisphere and :math:`-` the
downward one, and an operator with two superscripts maps its second
(incident) hemisphere to its first (outgoing) one.

Single-scattering properties
****************************

A layer has the particle extinction :math:`k_p` and scattering coefficient
:math:`\sigma` of its scattering set, and the scalar, unpolarized gas
extinction :math:`k_g\ge0`.  With the layer thickness :math:`\Delta z`,

.. math::

   k = k_p + k_g,\qquad
   \omega = \frac{\sigma}{k},\qquad
   \Delta\tau = k\,\Delta z.

The extinction is scalar: RT3 has no polarized extinction or dichroism.

The scattering-plane phase matrix of randomly oriented particles with a plane
of symmetry is that of :doc:`concept.tmatrix`,

.. math::

   \boldsymbol F(\Theta) =
   \begin{pmatrix}
   F_{11} & F_{12} & 0 & 0 \\
   F_{12} & F_{22} & 0 & 0 \\
   0 & 0 & F_{33} & F_{34} \\
   0 & 0 & -F_{34} & F_{44}
   \end{pmatrix},\qquad
   \frac{1}{4\pi}\int_{4\pi} F_{11}\,\mathrm d\Omega = 1,

six elements.  RT3 takes it as a Legendre series in :math:`\cos\Theta` of
such matrices,

.. math::

   \boldsymbol F(\Theta) = \sum_{l=0}^{L}\boldsymbol c_l\,P_l(\cos\Theta),
   \qquad (\boldsymbol c_0)_{11} = 1,

each coefficient again of six elements.  The coefficients include the factor
:math:`2l+1`: in DISORT's moments, :math:`(\boldsymbol c_l)_{11}=(2l+1)\chi_l`.
Rayleigh scattering is :math:`\boldsymbol c_0=(F_{11},F_{12},F_{22})=(1,-\tfrac12,1)`,
:math:`\boldsymbol c_1=(F_{33},F_{44})=(\tfrac32,\tfrac32)` and
:math:`\boldsymbol c_2=(F_{11},F_{12},F_{22})=(\tfrac12,\tfrac12,\tfrac12)`,
the other elements zero.

The phase matrix between the meridional bases of an incident direction
:math:`(\mu_i,\phi_i)` and a scattered direction :math:`(\mu_o,\phi_o)` is

.. math::

   \boldsymbol Z(\mu_o,\mu_i,u)
   = \boldsymbol L(\psi_o)\,\boldsymbol F(\Theta)\,\boldsymbol L(\psi_i),
   \qquad u = \phi_i - \phi_o,

   \cos\Theta = \mu_o\mu_i + \sqrt{(1-\mu_o^2)(1-\mu_i^2)}\cos u,

with the Stokes rotation

.. math::

   \boldsymbol L(\psi) =
   \begin{pmatrix}
   1 & 0 & 0 & 0\\
   0 & \cos2\psi & \sin2\psi & 0\\
   0 & -\sin2\psi & \cos2\psi & 0\\
   0 & 0 & 0 & 1
   \end{pmatrix}

:math:`\psi_i` and :math:`\psi_o` are the
angles between the scattering plane and the meridional planes of the two
directions.  Exactly forward and backward (:math:`\sin\Theta=0`) the
scattering plane is undefined and both rotations are the identity.  With
:math:`c_k=\cos2\psi_k` and :math:`s_k=\sin2\psi_k`, the product is

.. math::

   \boldsymbol Z =
   \begin{pmatrix}
   F_{11} & c_iF_{12} & s_iF_{12} & 0\\
   c_oF_{12} & c_ic_oF_{22}-s_is_oF_{33} & s_ic_oF_{22}+c_is_oF_{33} & s_oF_{34}\\
   -s_oF_{12} & -c_is_oF_{22}-s_ic_oF_{33} & -s_is_oF_{22}+c_ic_oF_{33} & c_oF_{34}\\
   0 & s_iF_{34} & -c_iF_{34} & F_{44}
   \end{pmatrix},

which is ARTS's laboratory-frame phase matrix of the same
:math:`\boldsymbol F` with the directions mapped as above.

Partition Stokes space into :math:`A=(I,Q)` and :math:`B=(U,V)` and let

.. math::

   \boldsymbol D = \operatorname{diag}(1,1,-1,-1),

the Mueller matrix of a mirror reflection, which changes the signs of
:math:`U` and :math:`V`.  Reversing the azimuth and
negating both direction cosines are both mirror reflections of the geometry:

.. math::

   \boldsymbol Z(\mu_o,\mu_i,-u)
   &= \boldsymbol D\,\boldsymbol Z(\mu_o,\mu_i,u)\,\boldsymbol D,\\
   \boldsymbol Z(-\mu_o,-\mu_i,u)
   &= \boldsymbol D\,\boldsymbol Z(\mu_o,\mu_i,u)\,\boldsymbol D.

So the diagonal :math:`2\times2` blocks :math:`\boldsymbol Z_{AA}` and
:math:`\boldsymbol Z_{BB}` are even in :math:`u` and the off-diagonal blocks
:math:`\boldsymbol Z_{AB}` and :math:`\boldsymbol Z_{BA}` odd.

Polarized radiative transfer
****************************

The equation solved in a layer is

.. math::

   \mu\frac{\partial\boldsymbol I}{\partial\tau}
   = \boldsymbol I
   - \frac{\omega}{4\pi}\int_{4\pi}
       \boldsymbol Z(\Omega,\Omega')\boldsymbol I(\Omega')
       \,\mathrm{d}\Omega'
   - \boldsymbol S_{\mathrm{dir}}
   - (1-\omega)B(\tau)\,\boldsymbol e_I,

with :math:`\boldsymbol e_I=[1,0,0,0]^{\mathsf T}`, the Planck function
:math:`B` of the layer temperature (linear in :math:`\tau` within a layer, see
`Sources and boundary conditions`_), and the direct-beam pseudo-source

.. math::

   \boldsymbol S_{\mathrm{dir}}(\tau,\Omega)
   = \frac{\omega}{4\pi}\,\frac{F_0}{\mu_0}\,e^{-\tau/\mu_0}\,
     \boldsymbol Z(\Omega,\Omega_0)\,\boldsymbol e_I,
   \qquad \Omega_0=(-\mu_0,0),

where :math:`F_0` is the direct flux on the horizontal at the top, so that
:math:`F_0/\mu_0` is the irradiance normal to the beam.

Azimuthal Fourier modes
=======================

The beam lies in the plane :math:`\phi=0` and the medium has the mirror
symmetry above.  The field is therefore symmetric under
:math:`\phi\to-\phi` with :math:`U` and :math:`V` odd, and RT3 expands

.. math::

   \boldsymbol I_A(\tau,\mu,\phi) = \sum_{m=0}^{M}\boldsymbol I^m_A(\tau,\mu)\cos m\phi,
   \qquad
   \boldsymbol I_B(\tau,\mu,\phi) = \sum_{m=0}^{M}\boldsymbol I^m_B(\tau,\mu)\sin m\phi,

with :math:`\boldsymbol I^0_B=0` and :math:`M` the azimuth order.  In the
terms of :doc:`concept.disort`, this is the one combined system that the beam
excites; the other has no source and vanishes.  Without the beam, only
:math:`m=0` has a source, and every :math:`m>0` mode and :math:`U` and
:math:`V` are zero.

Inserting the expansion, the scattering integral of mode :math:`m` couples
only mode :math:`m`.  With the azimuthal projections

.. math::

   \boldsymbol P^m_{AA} &= \frac{1}{2\pi}\int_0^{2\pi}\boldsymbol Z_{AA}(u)\cos mu\,\mathrm du,\\
   \boldsymbol P^m_{AB} &= \frac{1}{2\pi}\int_0^{2\pi}\boldsymbol Z_{AB}(u)\sin mu\,\mathrm du,\\
   \boldsymbol P^m_{BA} &= -\frac{1}{2\pi}\int_0^{2\pi}\boldsymbol Z_{BA}(u)\sin mu\,\mathrm du,\\
   \boldsymbol P^m_{BB} &= \frac{1}{2\pi}\int_0^{2\pi}\boldsymbol Z_{BB}(u)\cos mu\,\mathrm du,

(:math:`\boldsymbol P^0_{AB}=\boldsymbol P^0_{BA}=0` by the parity of the
blocks), the scattering term of mode :math:`m` is

.. math::

   \frac{\omega}{4\pi}\int_{4\pi}\boldsymbol Z\boldsymbol I\,\mathrm d\Omega'
   \;\to\;
   \frac{\omega}{2}\int_{-1}^{1}\boldsymbol P^m(\mu,\mu')\,
     \boldsymbol I^m(\mu')\,\mathrm d\mu'.

The direct beam gives the mode

.. math::

   \boldsymbol S^m_{\mathrm{dir}}(\tau,\mu)
   = \frac{\omega}{4\pi}\,\frac{F_0}{\mu_0}\,e^{-\tau/\mu_0}\,
     \boldsymbol z^m(\mu),

where :math:`\boldsymbol z^m` holds the azimuthal Fourier coefficients of the
first column :math:`\boldsymbol z(\mu,\phi)=\boldsymbol Z(\mu,\phi;-\mu_0,0)\boldsymbol e_I`:
the mean for :math:`m=0` (:math:`U` and :math:`V` zero), and for :math:`m>0`

.. math::

   \boldsymbol z^m_A = \frac1\pi\int_0^{2\pi}\boldsymbol z_A(\mu,\phi)\cos m\phi\,\mathrm d\phi,
   \qquad
   \boldsymbol z^m_B = \frac1\pi\int_0^{2\pi}\boldsymbol z_B(\mu,\phi)\sin m\phi\,\mathrm d\phi.

The mirror symmetries carry over to the modes,

.. math::

   \boldsymbol P^m(-\mu_o,-\mu_i) = \boldsymbol D\,\boldsymbol P^m(\mu_o,\mu_i)\,\boldsymbol D.

Discrete azimuth integrals
==========================

RT3 evaluates :math:`\boldsymbol Z` at :math:`N_\phi` azimuths
:math:`u_k=2\pi k/N_\phi`, of which it computes :math:`k\le N_\phi/2` and
takes the others from :math:`\boldsymbol Z(-u)=\boldsymbol D\boldsymbol Z(u)\boldsymbol D`,
and obtains the projections by a real FFT.  For a series of degree
:math:`L`,

.. math::

   N_\phi =
   \begin{cases}
     2\cdot2^{\lfloor\log_2(L+4)\rfloor+1}, & M>0,\\[2pt]
     2\lfloor(L+1)/2\rfloor+4, & M=0,
   \end{cases}

and for :math:`M=0` the mean is the plain average of the samples.
:math:`Z_{11}=F_{11}(\Theta)` is a polynomial of degree :math:`L` in
:math:`\cos u`, so its modes are sampled without aliasing.  The other elements
also depend on the rotation angles, which are not band-limited in :math:`u`;
their modes carry the aliasing of the sampling.

Stream equations
================

On the quadrature, with the signed streams :math:`\mu_i` of both hemispheres
and the repeated weights :math:`\overline w_j`, mode :math:`m` is the system

.. math::

   \mu_i\frac{\mathrm d\boldsymbol I^m_i}{\mathrm d\tau}
   = \boldsymbol I^m_i
   - \omega\sum_j \boldsymbol Z^m_{ij}\boldsymbol I^m_j
   - \boldsymbol S^m_{\mathrm{dir}}(\tau,\mu_i)
   - \delta_{m0}(1-\omega)B(\tau)\,\boldsymbol e_I,
   \qquad
   \boldsymbol Z^m_{ij} = \frac{\overline w_j}{2}\,\boldsymbol P^m(\mu_i,\mu_j),

in which a Mueller block maps the incident Stokes components (columns) to the
outgoing ones (rows).  Split by hemispheres, :math:`\boldsymbol Z^{m,++}`,
:math:`\boldsymbol Z^{m,+-}`, :math:`\boldsymbol Z^{m,-+}` and
:math:`\boldsymbol Z^{m,--}`, of which RT3 forms the first two and takes the
others from :math:`\boldsymbol Z^{m,--}=\boldsymbol D\boldsymbol Z^{m,++}\boldsymbol D`
and :math:`\boldsymbol Z^{m,-+}=\boldsymbol D\boldsymbol Z^{m,+-}\boldsymbol D`.

The discrete phase function must be normalized on the streams:

.. math::

   \frac12\sum_{i=1}^{N} w_i
   \left[\boldsymbol P^0(\mu_i,\mu_j)+\boldsymbol P^0(-\mu_i,\mu_j)\right]_{II} = 1

for every incident quadrature stream :math:`\pm\mu_j`, to :math:`10^{-7}` in
RT3.  The series must therefore not exceed the degree the quadrature
integrates,

.. math::

   L_{\max} =
   \begin{cases}
     4N-3, & \text{Gauss},\\
     2N-3, & \text{double Gauss},\\
     4N-5, & \text{Lobatto},
   \end{cases}
   \qquad L_{\max}\ge1,

to which RT3 truncates longer series; ARTS rejects a problem in which that
would drop a non-zero coefficient.

.. _Sec RT3 doubling-adding:

Doubling and adding
*******************

Instead of an eigen-decomposition, RT3 represents a slab by the operators that
map the radiance incident on it to the radiance leaving it.  For a slab with
the downward radiance :math:`\boldsymbol I^-_{\mathrm t}` incident at its top
and the upward radiance :math:`\boldsymbol I^+_{\mathrm b}` incident at its
bottom, the radiances leaving it are

.. math::

   \boldsymbol I^+_{\mathrm t}
   &= \boldsymbol R^{+-}\boldsymbol I^-_{\mathrm t}
    + \boldsymbol T^{++}\boldsymbol I^+_{\mathrm b} + \boldsymbol S^+,\\
   \boldsymbol I^-_{\mathrm b}
   &= \boldsymbol T^{--}\boldsymbol I^-_{\mathrm t}
    + \boldsymbol R^{-+}\boldsymbol I^+_{\mathrm b} + \boldsymbol S^-.

The operators are :math:`n\times n` matrices over the streams and Stokes
components, with :math:`n` the number of streams per hemisphere times
:math:`n_s`, and the sources are vectors of length :math:`n`, all per Fourier
mode.

Initial layer
=============

A scattering layer of optical thickness :math:`\Delta\tau` is halved
:math:`n_d` times to an initial sublayer of :math:`\delta\tau` at most
:math:`\delta\tau_{\max}`,

.. math::

   n_d = \max\!\left(0,\;
     \left\lfloor\log_2\frac{\max(\Delta\tau,10^{-7})}{\delta\tau_{\max}}\right\rfloor+1\right),
   \qquad
   \delta\tau = \Delta\tau\,2^{-n_d}.

:math:`\delta\tau` is a vertical optical thickness; along a stream the
sublayer is :math:`\delta\tau/\mu_i` thick.  To first order in
:math:`\delta\tau`, the sublayer transmits, reflects and emits

.. math::

   \boldsymbol T^{\pm\pm}_{ij}
   &= \delta_{ij}\boldsymbol 1
    - \frac{\delta\tau}{\mu_i}\left(\delta_{ij}\boldsymbol 1
    - \omega\,\boldsymbol Z^{m,\pm\pm}_{ij}\right),\\
   \boldsymbol R^{\pm\mp}_{ij}
   &= \frac{\delta\tau}{\mu_i}\,\omega\,\boldsymbol Z^{m,\pm\mp}_{ij},\\
   \boldsymbol S^{\pm}_i
   &= \frac{\delta\tau}{\mu_i}\,\boldsymbol J^m(\pm\mu_i),

with :math:`\mu_i>0` the cosine of the outgoing stream and
:math:`\boldsymbol J^m` the source of the sublayer: the direct-beam
pseudo-source at its top, or :math:`(1-\omega)B\,\boldsymbol e_I` in mode 0.

Adding
======

Slab 1 on top of slab 2 combine to a slab with

.. math::

   \boldsymbol\Gamma_+ &= \left(\boldsymbol 1 - \boldsymbol R_2^{+-}\boldsymbol R_1^{-+}\right)^{-1},\\
   \boldsymbol\Gamma_- &= \left(\boldsymbol 1 - \boldsymbol R_1^{-+}\boldsymbol R_2^{+-}\right)^{-1},\\
   \boldsymbol R^{+-} &= \boldsymbol R_1^{+-}
     + \boldsymbol T_1^{++}\boldsymbol\Gamma_+\boldsymbol R_2^{+-}\boldsymbol T_1^{--},\\
   \boldsymbol R^{-+} &= \boldsymbol R_2^{-+}
     + \boldsymbol T_2^{--}\boldsymbol\Gamma_-\boldsymbol R_1^{-+}\boldsymbol T_2^{++},\\
   \boldsymbol T^{++} &= \boldsymbol T_1^{++}\boldsymbol\Gamma_+\boldsymbol T_2^{++},\\
   \boldsymbol T^{--} &= \boldsymbol T_2^{--}\boldsymbol\Gamma_-\boldsymbol T_1^{--},\\
   \boldsymbol S^{+} &= \boldsymbol S_1^{+}
     + \boldsymbol T_1^{++}\boldsymbol\Gamma_+\left(\boldsymbol S_2^{+}+\boldsymbol R_2^{+-}\boldsymbol S_1^{-}\right),\\
   \boldsymbol S^{-} &= \boldsymbol S_2^{-}
     + \boldsymbol T_2^{--}\boldsymbol\Gamma_-\left(\boldsymbol S_1^{-}+\boldsymbol R_1^{-+}\boldsymbol S_2^{+}\right).

:math:`\boldsymbol\Gamma_+` and :math:`\boldsymbol\Gamma_-` sum the reflections
back and forth between the two slabs of the upward and the downward radiance
at their interface.

Doubling
========

Doubling adds a slab to an identical copy of itself, :math:`n_d` times, from
the initial sublayer to the layer.  The operators follow from the adding
formulas with both slabs equal.  The sources of the lower copy differ:

* the direct beam is attenuated across the upper copy, so its sources are
  :math:`e\boldsymbol S^\pm` with :math:`e=\exp(-\delta\tau'/\mu_0)` for the
  current slab thickness :math:`\delta\tau'`, and :math:`e\to e^2` at every
  doubling;
* the Planck function increases across the upper copy.  A second pair of
  vectors :math:`\boldsymbol C^\pm`, the source of the slab for a constant
  Planck function equal to that at its top (initially
  :math:`\boldsymbol C^\pm=\boldsymbol S^\pm`), gives the sources of the
  lower copy as :math:`\boldsymbol S^\pm+f\boldsymbol C^\pm`, where
  :math:`f` is the relative increase of the Planck function across the
  current slab, initially

  .. math::

     f = \frac{B_{\mathrm b}/B_{\mathrm t}-1}{2^{n_d}},

  for the Planck function :math:`B_{\mathrm t}` at the top and
  :math:`B_{\mathrm b}` at the bottom of the layer.  :math:`\boldsymbol C^\pm`
  is doubled as a source of its own, and :math:`f\to2f` at every doubling.

The thermal source of the layer is thereby that of a staircase: each initial
sublayer emits with the Planck function at its top, and the steps follow the
linear profile.

With :math:`\boldsymbol Z^{m,--}=\boldsymbol D\boldsymbol Z^{m,++}\boldsymbol D`
and :math:`\boldsymbol Z^{m,-+}=\boldsymbol D\boldsymbol Z^{m,+-}\boldsymbol D`,
every slab satisfies :math:`\boldsymbol T^{--}=\boldsymbol D\boldsymbol T^{++}\boldsymbol D`
and :math:`\boldsymbol R^{-+}=\boldsymbol D\boldsymbol R^{+-}\boldsymbol D`
(with :math:`\boldsymbol D` acting on every stream).  For :math:`n_s\le2`,
:math:`\boldsymbol D` is the identity on the transported components, the two
directions coincide, and RT3 computes them once.

.. _Sec RT3 absorbing layers:

Absorbing layers
================

A layer without scattering is integrated exactly.  With the slant optical
thickness :math:`p=\Delta\tau/\mu_i` and the Planck function :math:`B_{\mathrm t}`
at the top and :math:`B_{\mathrm b}` at the bottom, linear in :math:`\tau`,
:math:`s=(B_{\mathrm b}-B_{\mathrm t})/p`, in mode 0 and for :math:`I`,

.. math::

   \boldsymbol T^{\pm\pm} &= e^{-p}\,\boldsymbol 1,\qquad \boldsymbol R^{\pm\mp}=0,\\
   S^+_i &= B_{\mathrm t} + s - \left(B_{\mathrm t} + s(1+p)\right)e^{-p},\\
   S^-_i &= B_{\mathrm b} - s - \left(B_{\mathrm b} - s(1+p)\right)e^{-p}.

Sources and boundary conditions
*******************************

* The thermal source is unpolarized and enters mode 0 only: the Planck
  function (in W m\ :sup:`-2` Hz\ :sup:`-1` sr\ :sup:`-1`)
  of the temperatures at the layer boundaries, linear in :math:`\tau` within
  the layer (see `Doubling`_), times :math:`1-\omega`.
* The direct flux at level :math:`k` is
  :math:`F_0\exp(-\tau_k/\mu_0)`, with :math:`\tau_k` the total optical depth
  of the layers above it (delta-M scaled with delta-M).  Each initial sublayer
  takes the beam at its own top.
* The radiance incident at the top is isotropic, unpolarized blackbody
  radiation of the sky temperature, :math:`\boldsymbol I^-_{\mathrm{top}}=B(T_{\mathrm{sky}})\boldsymbol e_I`,
  in mode 0.

Surface reflection
******************

The surface is a slab of its own below the atmosphere, with
:math:`\boldsymbol R^{+-}=\boldsymbol R_{\mathrm s}`,
:math:`\boldsymbol T^{++}=\boldsymbol T^{--}=\boldsymbol 1`,
:math:`\boldsymbol R^{-+}=0` and no source, under which the surface radiance
:math:`\boldsymbol I_{\mathrm g}` is incident from below.  So at the surface

.. math::

   \boldsymbol I^+ = \boldsymbol R_{\mathrm s}\boldsymbol I^- + \boldsymbol I_{\mathrm g},

in every mode.

**Lambertian** (albedo :math:`A`), mode 0 only,

.. math::

   (\boldsymbol R_{\mathrm s})_{ij} = 2A\,\mu_jw_j\,\boldsymbol e_I\boldsymbol e_I^{\mathsf T},
   \qquad
   \boldsymbol I_{\mathrm g} =
     \left[(1-A)B(T_{\mathrm s}) + \frac{A}{\pi}F_0e^{-\tau_L/\mu_0}\right]\boldsymbol e_I,

the emission with the thermal source and the reflected beam with the solar
one.  It conserves energy on the streams when :math:`2\sum_j\mu_jw_j=1`, which
holds for the double-Gauss rule.

**Fresnel** (a plane interface to a medium of refractive index :math:`n`),
specular, in every mode,

.. math::

   (\boldsymbol R_{\mathrm s})_{ij} = \delta_{ij}\boldsymbol R_F(\mu_i),
   \qquad
   \boldsymbol R_F =
   \begin{pmatrix}
   R_1 & R_2 & 0 & 0\\
   R_2 & R_1 & 0 & 0\\
   0 & 0 & R_3 & -R_4\\
   0 & 0 & R_4 & R_3
   \end{pmatrix},

with the amplitude reflection coefficients :math:`r_v` and :math:`r_h`,
:math:`R_1=(|r_v|^2+|r_h|^2)/2`, :math:`R_2=(|r_v|^2-|r_h|^2)/2`,
:math:`R_3=\Re(r_vr_h^*)` and :math:`R_4=\Im(r_vr_h^*)`
and the emission, in mode 0,

.. math::

   \boldsymbol I_{\mathrm g} = (\boldsymbol 1-\boldsymbol R_F)B(T_{\mathrm s})\boldsymbol e_I
   = [(1-R_1)B,\,-R_2B,\,0,\,0]^{\mathsf T}.

The reflected direct beam would leave at :math:`\mu_0`, generally not a
stream, so a Fresnel surface is not allowed with the beam.

.. _Sec RT3 levels:

Radiances at the levels
***********************

For the radiance at a level, the slabs above it are added into one slab
:math:`a`, and those below it, with the surface, into one slab :math:`b`.
The radiances at the level then follow from the two slabs and the radiances
incident on the whole,

.. math::

   \boldsymbol I^- &= \left(\boldsymbol 1-\boldsymbol R_a^{-+}\boldsymbol R_b^{+-}\right)^{-1}
     \left[\boldsymbol T_a^{--}\boldsymbol I^-_{\mathrm{top}} + \boldsymbol S_a^-
     + \boldsymbol R_a^{-+}\left(\boldsymbol T_b^{++}\boldsymbol I_{\mathrm g}+\boldsymbol S_b^+\right)\right],\\
   \boldsymbol I^+ &= \left(\boldsymbol 1-\boldsymbol R_b^{+-}\boldsymbol R_a^{-+}\right)^{-1}
     \left[\boldsymbol T_b^{++}\boldsymbol I_{\mathrm g} + \boldsymbol S_b^+
     + \boldsymbol R_b^{+-}\left(\boldsymbol T_a^{--}\boldsymbol I^-_{\mathrm{top}}+\boldsymbol S_a^-\right)\right].

The radiance at the azimuth :math:`\phi` is the sum of the modes,
:math:`\sum_m\boldsymbol I^m_A\cos m\phi` for :math:`I` and :math:`Q` and
:math:`\sum_m\boldsymbol I^m_B\sin m\phi` for :math:`U` and :math:`V`.

Delta-M scaling
***************

With delta-M, every scattering set is scaled before the transport with the
classical delta-M of :doc:`concept.disort` (:math:`r_l=1`), applied to the
whole matrix.  The truncation order is :math:`2N_{\mathrm{tot}}`, with
:math:`N_{\mathrm{tot}}` the number of streams per hemisphere including the
extra directions, and

.. math::

   f &= \frac{(\boldsymbol c_{2N_{\mathrm{tot}}})_{11}}{4N_{\mathrm{tot}}+1},\qquad
   s = 1-\omega f,\\
   k' &= s\,k,\qquad
   \omega' = \frac{\omega(1-f)}{1-\omega f},\\
   \boldsymbol c'_l &= (2l+1)\,\frac{\boldsymbol c_l/(2l+1) - f\,\boldsymbol 1}{1-f},
   \qquad l<2N_{\mathrm{tot}},

where :math:`\boldsymbol 1` is the identity in compact form
(:math:`F_{11}=F_{22}=F_{33}=F_{44}=1`): the forward peak removed is
:math:`f` times the identity, so the diagonal elements lose :math:`f` and
:math:`F_{12}` and :math:`F_{34}` are only renormalized.  The series is then truncated to
:math:`L_{\max}`.  The scaled extinction attenuates the beam; the thermal
emission of a layer is unchanged, as
:math:`(1-\omega')\,\Delta\tau'=(1-\omega)\,\Delta\tau`.  RT3 applies no
correction for the removed peak (no IMS or TMS).

Arbitrary-angle radiances
*************************

Directions other than the quadrature nodes are added as streams of weight
zero.  They receive scattering from the quadrature streams and are reflected
by the surface, but scatter nothing back, so the radiances on the nodes do
not change.  Their radiances are the solution along those directions of the
same discrete source function, to the accuracy of the doubling.

Fluxes
******

Only mode 0 contributes to the hemispheric fluxes.  RT3 returns them for
every Stokes component,

.. math::

   \boldsymbol F^\uparrow = 2\pi\sum_{j=1}^N w_j\mu_j \boldsymbol I^{0,+}_j,\qquad
   \boldsymbol F^\downarrow = 2\pi\sum_{j=1}^N w_j\mu_j \boldsymbol I^{0,-}_j
     + F_0e^{-\tau/\mu_0}\boldsymbol e_I,

with the direct flux added to the downward :math:`I`.

Mathematical limitations
************************

* **First-order initial layer.**  The thin-layer operators are first order in
  :math:`\delta\tau`, so the error of a scattering layer is first order in
  :math:`\delta\tau_{\max}`, and larger at small :math:`\mu_i`, where the
  slant thickness is larger.  Absorbing layers are exact.
* **Thermal staircase.**  Within each initial sublayer the Planck function
  is that at its top, an error of first order in :math:`\delta\tau`.
* **Azimuth sampling.**  The modes of the polarized elements alias at the
  sampling of `Discrete azimuth integrals`_.
* **Legendre degree.**  The series is limited to :math:`L_{\max}` of the
  quadrature; a sharply peaked phase function needs delta-M or more streams.
* **Plane-parallel.**  There is no pseudo-spherical beam.
* **Surfaces.**  A Fresnel surface cannot reflect the beam, and a Lambertian
  surface conserves energy on the streams only for the double-Gauss rule.

See :doc:`dev.rt3` for the implementation and its validation.
