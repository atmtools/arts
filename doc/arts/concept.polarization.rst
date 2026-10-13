.. _Sec Polarization geometry:

Polarization geometry
#####################

A Stokes vector describes a field only relative to a pair of reference
directions across its ray.  This page defines the pair that ARTS uses, how a
Stokes vector changes when the pair changes, and how scattering and surface
reflection are built from it.  The Stokes components themselves are
introduced in :ref:`the Stokes vector section <stok-vec>`.

.. _Sec Polarization basis:

The polarization basis of a ray
*******************************

A ray is described by its propagation direction :math:`\hat k`, the direction
in which the radiation travels.  A ray that goes down has a downward
:math:`\hat k`, also when the calculation follows it backwards from a
sensor.  With :math:`\hat z` the local up direction, the meridional basis of
the ray is

.. math::

   \hat h = \frac{\hat k\times\hat z}{|\hat k\times\hat z|},\qquad
   \hat v = \hat h\times\hat k .

:math:`\hat h` is horizontal, :math:`\hat v` lies in the vertical plane of
:math:`\hat k` on its upward side (:math:`\hat v\cdot\hat z\geq0`), and
:math:`(\hat v,\hat h,\hat k)` is right-handed.  The same rule holds in both
hemispheres.  In the local east-north-up frame, a propagation direction with
zenith angle :math:`\theta` and azimuth :math:`\alpha`, clockwise from north
seen from above, is

.. math::

   \hat k = (\sin\theta\sin\alpha,\ \sin\theta\cos\alpha,\ \cos\theta),
   \qquad
   \hat h = (\cos\alpha,\ -\sin\alpha,\ 0).

For a ray along :math:`\pm\hat z` the cross product vanishes.  Its basis is
then the limit along its meridian, the same :math:`\hat h` as above, so it is
fixed by the ray's azimuth.

With the field :math:`\vec E = E_v\hat v + E_h\hat h`, the Stokes vector in
this basis is

.. math::

   I = |E_v|^2+|E_h|^2,\qquad
   Q = |E_v|^2-|E_h|^2,\qquad
   U = 2\,\Re\left(E_vE_h^*\right),\qquad
   V = 2\,\Im\left(E_vE_h^*\right).

Positive :math:`Q` is vertical polarization.  With the time dependence
:math:`\vec E(t) = \Re\left[(E_v\hat v + E_h\hat h)e^{-i\omega t}\right]`,
positive :math:`V` is a field that rotates counterclockwise seen looking in
the direction of propagation: left circular polarization, as defined in
:ref:`the Stokes vector section <stok-vec>`.  These definitions are those of
the polarized discrete-ordinate solvers (:doc:`concept.disort`,
:doc:`concept.rt3`, :doc:`concept.rt4`) and of the scattering and reflection
geometry below.

Changing the basis
******************

Two bases of the same ray that are both right-handed with :math:`\hat k`
differ by a rotation about :math:`\hat k` through the angle :math:`\psi`, with

.. math::

   \cos\psi = \hat v'\cdot\hat v,\qquad \sin\psi = \hat v'\cdot\hat h .

The Stokes vector in the new basis is :math:`\vec I' = \mathbf L(\psi)\vec I`,

.. math::

   \mathbf L(\psi) =
   \begin{bmatrix}
      1 & 0 & 0 & 0\\
      0 & \cos2\psi & \sin2\psi & 0\\
      0 & -\sin2\psi & \cos2\psi & 0\\
      0 & 0 & 0 & 1
   \end{bmatrix}.

:math:`I` and :math:`V` do not change.  A basis of the other handedness, such
as :math:`(\hat v,-\hat h)`, is not reached by a rotation: it changes the
signs of :math:`U` and :math:`V` together.  The two handednesses are never
mixed.

Near the zenith and the nadir, a small change of direction can turn the
meridional basis by a large angle.  Two rays there that differ by a small
angle but have different azimuths have bases rotated by the difference of
their azimuths.

Jones and Mueller matrices
**************************

A linear interaction that does not depolarize maps the field components of
the incident ray to those of the outgoing ray, each in its own basis, by a
complex :math:`2\times2` Jones matrix,

.. math::

   \begin{bmatrix}E_v'\\E_h'\end{bmatrix}
   = \mathbf J\begin{bmatrix}E_v\\E_h\end{bmatrix}.

Its Mueller matrix follows from the Stokes definitions above,

.. math::

   \mathbf M = \mathbf A\left(\mathbf J\otimes\mathbf J^*\right)\mathbf A^{-1},
   \qquad
   \mathbf A =
   \begin{bmatrix}
      1 & 0 & 0 & 1\\
      1 & 0 & 0 & -1\\
      0 & 1 & 1 & 0\\
      0 & -i & i & 0
   \end{bmatrix},

where :math:`\mathbf A` takes the products
:math:`(E_vE_v^*, E_vE_h^*, E_hE_v^*, E_hE_h^*)` to the Stokes vector.  Its
determinant, :math:`|\det\mathbf J|^4`, is not negative, and changes of basis
keep it so.  A Mueller matrix with a negative determinant therefore describes
no non-depolarizing interaction in any choice of bases.

.. _Sec Polarization scattering:

Scattering
**********

For randomly oriented particles with a plane of symmetry, the scattering
matrix :math:`\mathbf F(\Theta)` is defined in the scattering-plane bases of
the two rays :cite:p:`Mishchenko:02`.  With

.. math::

   \hat N = \frac{\hat k_\mathrm{in}\times\hat k_\mathrm{out}}
                 {|\hat k_\mathrm{in}\times\hat k_\mathrm{out}|},

the scattering-plane basis of a ray along :math:`\hat k` is
:math:`(\hat N\times\hat k, \hat N)`: its first vector lies in the scattering
plane and its second is perpendicular to it.  :math:`\mathbf F` has the block
form

.. math::

   \mathbf F =
   \begin{bmatrix}
      F_{11} & F_{12} & 0 & 0\\
      F_{12} & F_{22} & 0 & 0\\
      0 & 0 & F_{33} & F_{34}\\
      0 & 0 & -F_{34} & F_{44}
   \end{bmatrix}.

The phase matrix between the meridional bases of the two rays is

.. math::

   \mathbf Z = \mathbf L(\psi_\mathrm{out})\,\mathbf F(\Theta)\,
               \mathbf L(\psi_\mathrm{in}),

where :math:`\psi_\mathrm{in}` turns the incident meridional basis into the
incident scattering-plane basis and :math:`\psi_\mathrm{out}` the outgoing
scattering-plane basis into the outgoing meridional basis, each as a change
of basis above.

Exactly forward or backward, :math:`\hat k_\mathrm{in}\times\hat
k_\mathrm{out}=0` and the scattering plane is undefined.  The plane through
:math:`\hat k_\mathrm{in}` and its :math:`\hat h` is then used.  For a regular
scattering matrix, with :math:`F_{12}=F_{34}=0` and :math:`F_{22}=\pm F_{33}`
at :math:`\Theta=0` and :math:`\pi`, :math:`\mathbf Z` does not depend on that
choice.  The rotations are taken from the bases themselves also for nearly
forward and backward pairs, since near the zenith and the nadir the bases of
nearly parallel rays can differ by a large rotation.

The phase matrices on the streams of the polarized solvers, and the
laboratory-frame phase matrices of ARTS's randomly oriented scattering data,
are built this way.

.. _Sec Polarization reflection:

Surface reflection
******************

A surface element whose normal :math:`\hat n` bisects
:math:`-\hat k_\mathrm{in}` and :math:`\hat k_\mathrm{out}` reflects the one
ray into the other: a flat surface for the specular direction, or a facet of
a rough surface.  The reflection is described in the plane of incidence.

.. math::

   \hat s = \frac{\hat k_\mathrm{in}\times\hat n}
                 {|\hat k_\mathrm{in}\times\hat n|}

is perpendicular to that plane, and each wave has :math:`\hat p = \hat
s\times\hat k`, so :math:`(\hat p,\hat s,\hat k)` is right-handed for both.
The surface reflects the two components independently,

.. math::

   E_p' = r_v E_p,\qquad E_s' = r_h E_s,

with the Fresnel amplitudes

.. math::

   r_v = \frac{n_2\cos\theta_1-n_1\cos\theta_2}{n_2\cos\theta_1+n_1\cos\theta_2},
   \qquad
   r_h = \frac{n_1\cos\theta_1-n_2\cos\theta_2}{n_1\cos\theta_1+n_2\cos\theta_2},

where :math:`\theta_1` is the local angle of incidence, between
:math:`-\hat k_\mathrm{in}` and :math:`\hat n`, and :math:`\theta_2` the angle
of refraction.  A perfect conductor has :math:`r_v=1` and :math:`r_h=-1`.

The Mueller matrix of the reflection is the rotation from the incident
meridional basis to :math:`(\hat p,\hat s)`, the Fresnel matrix

.. math::

   \mathbf R_F = \frac12
   \begin{bmatrix}
      |r_v|^2+|r_h|^2 & |r_v|^2-|r_h|^2 & 0 & 0\\
      |r_v|^2-|r_h|^2 & |r_v|^2+|r_h|^2 & 0 & 0\\
      0 & 0 & 2\,\Re(r_vr_h^*) & -2\,\Im(r_vr_h^*)\\
      0 & 0 & 2\,\Im(r_vr_h^*) & 2\,\Re(r_vr_h^*)
   \end{bmatrix},

and the rotation from :math:`(\hat p,\hat s)` to the outgoing meridional
basis.  For a horizontal surface both rotations vanish, and the reflection is
:math:`\mathbf R_F` itself.

At normal incidence, :math:`\mathbf R_F` of a perfect conductor is
:math:`\operatorname{diag}(1,1,-1,-1)`.  The reversal of circular
polarization by a mirror is part of the reflection matrix, and the reflected
Stokes vector receives no further change of sign.  The surface emits what
makes up the difference to a blackbody, :math:`(\mathbf 1-\mathbf R)\vec B`
for the Planck vector :math:`\vec B`.

.. _Sec Polarization reflectance without phase:

Reflectances without phase
**************************

Some surface descriptions give only the power reflectances
:math:`|r_v|^2` and :math:`|r_h|^2`: a single scalar reflectance, or the
emissivities :math:`1-|r_v|^2` and :math:`1-|r_h|^2` of empirical models
such as TESSEM and TELSEM.  They lack the phase :math:`\delta` of

.. math::

   r_vr_h^* = |r_v||r_h|\,e^{i\delta},

which sets the :math:`U` and :math:`V` block of :math:`\mathbf R_F`: its
diagonal is :math:`|r_v||r_h|\cos\delta` and its coupling of :math:`U` and
:math:`V` is :math:`\mp|r_v||r_h|\sin\delta`.

At normal incidence a surface cannot tell :math:`\hat v` from :math:`\hat h`,
so :math:`r_v=-r_h` in the convention above and :math:`\delta=\pi`.  A
loss-free dielectric keeps :math:`\delta=\pi` up to the Brewster angle.  ARTS
therefore takes :math:`\delta=\pi` for these descriptions: the :math:`U` and
:math:`V` diagonal elements are :math:`-|r_v||r_h|` and the two are not
coupled.  A scalar reflectance :math:`R` is then the isotropic mirror
:math:`R\operatorname{diag}(1,1,-1,-1)`.  The intermediate phases of absorbing
surfaces, and the coupling of :math:`U` and :math:`V` that they cause, cannot
be represented by such descriptions.
