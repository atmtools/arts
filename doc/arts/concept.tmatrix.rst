T-matrix particle scattering
============================

The T-matrix method relates an incident electromagnetic field to the field
scattered by a particle.  Both fields are expanded in vector spherical waves;
for expansion-coefficient vectors :math:`\boldsymbol a` and
:math:`\boldsymbol b`, the linear relation is

.. math::

   \boldsymbol b = \boldsymbol T\boldsymbol a.

For a given particle and wavelength, the same T-matrix can be used for different
illumination and observation directions.  The particle model considered here
is a homogeneous spheroid, finite circular cylinder, Chebyshev particle or
generalized Chebyshev particle.

Particle size and orientation
-----------------------------

Volume-equivalent radius is the radius of a sphere with the same volume as the
particle.  Surface-area-equivalent radius is the radius of a sphere with the
same surface area.  These definitions agree for a sphere, but generally differ
for nonspherical particles.

A fixed-orientation calculation specifies the particle orientation as well as
the incident and scattered directions.  A random-orientation calculation
averages over particle orientations.  An azimuthally random calculation
averages over the azimuths of the particle's symmetry axis at given tilts from
the vertical, such as of falling particles that the flow aligns.  Averaging
over a particle size distribution is an additional operation.

Amplitude and phase matrices
----------------------------

The complex 2 by 2 Jones amplitude matrix acts on electric-field polarization
components.  The corresponding real 4 by 4 Mueller phase matrix acts on the
Stokes vector.  Its elements are quadratic combinations of the complex
amplitudes.  With an amplitude measured in length units, the phase matrix has
units of length squared.

For randomly oriented particles of the shapes considered here, the normalized
phase matrix in the scattering-plane basis has the form

.. math::

   \boldsymbol F(\Theta) =
   \begin{pmatrix}
   F_{11} & F_{12} & 0 & 0 \\
   F_{12} & F_{22} & 0 & 0 \\
   0 & 0 & F_{33} & F_{34} \\
   0 & 0 & -F_{34} & F_{44}
   \end{pmatrix},\qquad
   \int_{4\pi} F_{11}\,d\Omega = 4\pi.

These matrix elements are dimensionless.  Multiplication by
:math:`C_{\mathrm{sca}}/(4\pi)` gives the differential scattering cross-section
matrix, where :math:`C_{\mathrm{sca}}` is the scattering cross section.
For a size distribution, the cross sections and normalized phase matrix must
refer to the same distribution average.

The results are in ARTS's Stokes convention (:doc:`concept.polarization`).  The
equations below are those of :cite:t:`Mishchenko:00` and :cite:t:`Mishchenko:91`,
whose Stokes parameters, like those of Bohren and Huffman, have :math:`U` of
the other sign and the same :math:`I`, :math:`Q` and :math:`V`.  The amplitude
matrix needs no change; the phase matrix is its Mueller matrix in ARTS's
convention, and :math:`F_{34}` of random orientation is the negative of
Mishchenko's.

See :doc:`user.tmatrix` for ARTS parameter names, units, and result types.

.. _Sec T-matrix equations:

The equations solved
--------------------

This section states what the ARTS implementation computes, in the order it
computes it.  The method is the extended boundary condition method of
:cite:t:`Waterman:71` in the formulation of :cite:t:`Mishchenko:98`; the
fixed-orientation amplitude matrix follows :cite:t:`Mishchenko:00` and the
analytical orientation average :cite:t:`Mishchenko:91`.  Throughout,
:math:`\lambda` is the wavelength and :math:`k=2\pi/\lambda` the wavenumber in
the surrounding medium, and :math:`\eta=\eta_r+i\eta_i` is the refractive index
of the particle relative to that medium.  :math:`n` and :math:`n'` are
multipole degrees and :math:`m` an azimuthal order; :math:`\eta` is never an
index.

.. _Sec T-matrix shapes:

Particle shapes
^^^^^^^^^^^^^^^

The particle is rotationally symmetric about the :math:`z` axis of its own
frame.  Its surface is :math:`r(\theta)` in the polar angle :math:`\theta`, and

.. math::

   \rho(\theta) = \frac{1}{r(\theta)}\frac{dr(\theta)}{d\theta}.

The shapes are parameterized by the volume-equivalent radius :math:`r_v`.  A
**spheroid** with horizontal-to-rotational axis ratio :math:`\varepsilon` has
the horizontal semi-axis :math:`a=r_v\varepsilon^{1/3}` and

.. math::

   r^2(\theta) = \frac{a^2}{\sin^2\theta+\varepsilon^2\cos^2\theta},\qquad
   \rho(\theta) = \frac{(\varepsilon^2-1)\sin\theta\cos\theta}{\sin^2\theta+\varepsilon^2\cos^2\theta}.

A **cylinder** with diameter-to-length ratio :math:`\varepsilon` has the
half-length :math:`h=r_v\,(2/(3\varepsilon^2))^{1/3}` and the radius
:math:`\varepsilon h`.  Its edge is at :math:`\tan\theta_e=\varepsilon`, and for
:math:`0\le\theta\le\pi/2`

.. math::

   r(\theta) = \frac{h}{\cos\theta},\quad \rho(\theta)=\tan\theta
   \quad(\theta<\theta_e);\qquad
   r(\theta) = \frac{\varepsilon h}{\sin\theta},\quad \rho(\theta)=-\cot\theta
   \quad(\theta>\theta_e).

A **Chebyshev particle** :cite:p:`Mishchenko:98` of degree :math:`n\ge1` and
deformation :math:`\varepsilon`, :math:`|\varepsilon|<1`, has

.. math::

   r(\theta) = r_0\,(1+\varepsilon\cos n\theta),\qquad
   \rho(\theta) = -\frac{\varepsilon n\sin n\theta}{1+\varepsilon\cos n\theta},\qquad
   r_0 = r_v A^{-1/3},

with :math:`A=1+\tfrac32\varepsilon^2(4n^2-2)/(4n^2-1)` for odd :math:`n`, and
for even :math:`n` additionally
:math:`-3\varepsilon(1+\varepsilon^2/4)/(n^2-1)-\varepsilon^3/(4(9n^2-1))`.

A **generalized Chebyshev particle**, the original program's shape of
distorted water drops, has the coefficients :math:`c_0,\ldots,c_N` and

.. math::

   r(\theta) = r_0\,a(\theta),\quad a(\theta) = 1+c_0+\sum_{n=1}^{N}c_n\cos n\theta,\qquad
   \rho(\theta) = -\frac{1}{a(\theta)}\sum_{n=1}^{N}n\,c_n\sin n\theta,\qquad
   r_0 = r_v V^{-1/3},

with :math:`V=\tfrac12\int_{-1}^{1}a^3\,d\cos\theta`, the volume of
:math:`r=a` over that of the unit sphere.  The coefficients are the user's,
such as of the shape of a raindrop of a given size, where the original fixes
those of one drop; the Chebyshev particle is the case of the single
coefficient :math:`c_n=\varepsilon`.
:math:`\cos n\theta=T_n(\cos\theta)`, so :math:`a^3` is a polynomial of
degree :math:`3N` in :math:`\cos\theta`, and the
:math:`\lceil(3N+1)/2\rceil`-point Gauss-Legendre rule gives :math:`V`
exactly.  The coefficients must satisfy :math:`1+c_0>\sum_{n\ge1}|c_n|`,
which keeps :math:`r` positive.

Spheroids, cylinders, Chebyshev particles of even degree and generalized
Chebyshev particles without odd coefficients are also mirror symmetric in the
plane :math:`z=0`, :math:`r(\pi-\theta)=r(\theta)`, since
:math:`\cos n(\pi-\theta)=(-1)^n\cos n\theta`.  The others are not.

When the given radius is the surface-area-equivalent radius :math:`r_s`, it is
converted by :math:`r_v=q\,r_s` with

.. math::

   q &= \left[\tfrac12\left(\varepsilon^{2/3}+\varepsilon^{-1/3}\frac{\arcsin e}{e}\right)\right]^{-1/2},
   \quad e=\sqrt{1-\varepsilon^2} &&\text{(spheroid, } \varepsilon<1),\\
   q &= \left[\tfrac14\left(2\varepsilon^{2/3}+\varepsilon^{-4/3}\frac{1}{e}\ln\frac{1+e}{1-e}\right)\right]^{-1/2},
   \quad e=\sqrt{1-\varepsilon^{-2}} &&\text{(spheroid, } \varepsilon>1),\\
   q &= \left(\frac{3}{2\varepsilon}\right)^{1/3}\left(\frac{\varepsilon+2}{2\varepsilon}\right)^{-1/2}
   &&\text{(cylinder)}.

For :math:`\varepsilon=1` both spheroid expressions are :math:`0/0` with the
limit :math:`q=1`, and :math:`\ln\frac{1+e}{1-e}` is evaluated as
:math:`2\operatorname{artanh}e`, which does not cancel as :math:`e\to0`.

For a Chebyshev particle, with :math:`a=1+\varepsilon\cos n\theta` and
:math:`b=\varepsilon n\sin n\theta`, the volume of :math:`r=a` over that of
the unit sphere is :math:`A` above, exactly, and

.. math::

   q = \frac{A^{1/3}}{\bigl(\tfrac12\int_{-1}^{1}a\sqrt{a^2+b^2}\,d\cos\theta\bigr)^{1/2}},

the volume- over the surface-area-equivalent radius of :math:`r=a`.  A
generalized Chebyshev particle has the same with its :math:`a`,
:math:`b=-da/d\theta` and :math:`V` in place of :math:`A`.  The area
integral is not a polynomial in :math:`\cos\theta`.  It is integrated by
Boost's adaptive 61-point Gauss-Kronrod rule until the error estimate is at most
:math:`\delta'` relative to the integral, the criterion of
:ref:`Sec T-matrix convergence`; the conversion fails if the estimate does not
reach it within Boost's default depth of 15 bisections.  A sphere has
:math:`q=1`.

The maximum dimension is the largest distance between two points of the
surface.  Two points at opposite azimuths and polar angles :math:`\theta_1`,
:math:`\theta_2` are

.. math::

   D(\theta_1,\theta_2) = \sqrt{(r_1\sin\theta_1+r_2\sin\theta_2)^2+(r_1\cos\theta_1-r_2\cos\theta_2)^2}

apart, which no other azimuths exceed.  It is :math:`2r_v\max(\varepsilon^{1/3},\varepsilon^{-2/3})`
for a spheroid and the diagonal :math:`2h\sqrt{1+\varepsilon^2}` for a
cylinder.  For a Chebyshev particle of degree :math:`n`, and a generalized one
of :math:`N=n`, it is found on the grid
of :math:`4(n+1)+1` equidistant polar angles in each of :math:`\theta_1` and
:math:`\theta_2`: :math:`D` oscillates at most :math:`n+1` times over
:math:`[0,\pi]` in each, so four steps per half-oscillation put grid points on
both sides of every maximum.  Every local maximum of the grid is maximized
further by Brent's method in each angle in turn, within one grid step, until
:math:`D` no longer grows.

.. _Sec T-matrix quadrature:

Surface quadrature
^^^^^^^^^^^^^^^^^^

Surface integrals of a mirror-symmetric particle are taken over
:math:`x=\cos\theta` on one hemisphere, :math:`0<x<1`, with nodes :math:`x_i`
and weights :math:`w_i`, :math:`i=1,\ldots,N_G`; the mirror symmetry supplies
the other hemisphere (:ref:`Sec T-matrix integrals`).  For the
mirror-symmetric spheroid and (generalized) Chebyshev particles these are the
positive half of the :math:`2N_G`-point Gauss-Legendre rule.  The other
(generalized) Chebyshev particles use the whole :math:`2N_G`-point rule on
:math:`-1<x<1`.  For a cylinder the
:math:`\lfloor N_G/2\rfloor`-point Gauss-Legendre rule on
:math:`[\cos\theta_e,1]` (the end cap) and the
:math:`(N_G-\lfloor N_G/2\rfloor)`-point rule on :math:`[0,\cos\theta_e]` (the
side) together keep the kink at the edge out of the integrands.  At each node,

.. math::

   s_i=\sin\theta_i,\quad r_i=r(\theta_i),\quad \rho_i=\rho(\theta_i),\quad
   W_i = w_i r_i^2,\quad \xi_i = k r_i,\quad \zeta_i=\eta\,\xi_i.

.. _Sec T-matrix functions:

Special functions
^^^^^^^^^^^^^^^^^

:math:`j_n` and :math:`y_n` are the spherical Bessel functions of the first and
second kind, and :math:`h_n=j_n+iy_n`.  For any of them,

.. math::

   \tilde f_n(z) = \frac{1}{z}\frac{d}{dz}\bigl[z f_n(z)\bigr] = f_{n-1}(z)-\frac{n}{z}f_n(z).

At each node the functions are :math:`j_n(\xi_i)`, :math:`y_n(\xi_i)` and
:math:`j_n(\zeta_i)` for :math:`n=1,\ldots,n_{\max}`, with their tilde forms.
:math:`y_n` comes from upward recurrence.  Both :math:`j_n` come from downward
recurrence of the ratio :math:`j_n/j_{n-1}`, started at
:math:`j_L/j_{L-1}=z/(2L+1)` with :math:`L=n_{\max}+\Delta n`, and normalized
by :math:`j_{-1}(z)=\cos z/z`, that is :math:`j_0=(j_0/j_{-1})\cos z/z`, with

.. math::

   \Delta n_\xi = \bigl\lfloor 1.2\sqrt{\max(\xi_{\max},n_{\max})}\bigr\rfloor+3,\qquad
   \Delta n_\zeta = \bigl\lfloor t+4t^{1/3}+1.2\sqrt{t}\bigr\rfloor-n_{\max}+5,\quad
   t=\max(|\eta|\,\xi_{\max},n_{\max}),

where :math:`\xi_{\max}=\max_i\xi_i`.  These are the starting rules of the
double-precision original; the factor 1.2 is what extended precision would
have to raise.

The angular functions are the Wigner functions :math:`d^n_{0m}(\theta)`,
:math:`n\ge m\ge 0`.  :math:`d^n_{00}=P_n(\cos\theta)`, and for :math:`m>0`

.. math::

   d^n_{0m}(\theta) = \sqrt{\frac{(n-m)!}{(n+m)!}}\,P_n^m(\cos\theta)
   = \frac{1}{\sqrt 2}S_n^m(\cos\theta),

with :math:`P_n^m` without the Condon-Shortley phase and :math:`S_n^m` the
Schmidt semi-normalized functions.  With them,

.. math::

   \tau^m_n(\theta) = \frac{d}{d\theta}d^n_{0m}(\theta),\qquad
   \pi^m_n(\theta) = \frac{m\,d^n_{0m}(\theta)}{\sin\theta}.

At the poles, where :math:`\pi^m_n` is a limit, :math:`\pi^1_n` and
:math:`\tau^1_n` are :math:`\tfrac12\sqrt{n(n+1)}` at :math:`\theta=0`,
:math:`\pi^1_n(\pi)=(-1)^{n+1}\tfrac12\sqrt{n(n+1)}` and
:math:`\tau^1_n(\pi)=(-1)^n\tfrac12\sqrt{n(n+1)}`, and all other orders vanish.
The Schmidt recurrences behind :math:`d^n_{0m}` carry the factor
:math:`\sin^m\theta` and never divide by :math:`\sin\theta`, so with their
:math:`\sin\theta` the expression for :math:`\pi^m_n` is accurate to rounding
up to the poles.  It is :math:`0/0` only where :math:`\sin\theta=0`, and only
there does the amplitude matrix use the limits.

.. _Sec T-matrix integrals:

The T-matrix of one azimuthal order
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The T-matrix of a rotationally symmetric particle separates into one matrix per
azimuthal order :math:`m`.  For :math:`m\ge0` its degrees are
:math:`n,n'=n_0,\ldots,n_{\max}`, :math:`n_0=\max(m,1)`.  With
:math:`d_n=d^n_{0m}(\theta_i)`, :math:`d'_n=\tau^m_n(\theta_i)`,
:math:`h_n=h_n(\xi_i)`, :math:`J_{n'}=j_{n'}(\zeta_i)` and

.. math::

   a_{nn'} = \frac12\sqrt{\frac{(2n+1)(2n'+1)}{n(n+1)\,n'(n'+1)}},

the surface integrals are, for :math:`n+n'` odd,

.. math::

   K^{11}_{nn'} &= 2a_{nn'}\sum_i W_i\frac{m}{s_i}\,(d_nd'_{n'}+d'_nd_{n'})\,J_{n'}h_n,\\
   K^{22}_{nn'} &= 2a_{nn'}\sum_i W_i\frac{m}{s_i}\Bigl[(d_nd'_{n'}+d'_nd_{n'})\,\tilde J_{n'}\tilde h_n
     +\rho_i\,d_nd_{n'}\Bigl(\frac{n(n+1)}{\xi_i}\tilde J_{n'}h_n+\frac{n'(n'+1)}{\zeta_i}J_{n'}\tilde h_n\Bigr)\Bigr],

and for :math:`n+n'` even,

.. math::

   K^{12}_{nn'} &= 2a_{nn'}\sum_i W_i\Bigl[\Bigl(\frac{m^2}{s_i^2}d_nd_{n'}+d'_nd'_{n'}\Bigr)J_{n'}\tilde h_n
     +\rho_i\,\frac{n(n+1)}{\xi_i}\,d_nd'_{n'}\,J_{n'}h_n\Bigr],\\
   K^{21}_{nn'} &= 2a_{nn'}\sum_i W_i\Bigl[\Bigl(\frac{m^2}{s_i^2}d_nd_{n'}+d'_nd'_{n'}\Bigr)\tilde J_{n'}h_n
     +\rho_i\,\frac{n'(n'+1)}{\zeta_i}\,d'_nd_{n'}\,J_{n'}h_n\Bigr].

For a mirror-symmetric particle the other parity of each vanishes: under
:math:`\theta\to\pi-\theta` the integrand of :math:`K^{12}` and :math:`K^{21}`
changes by :math:`(-1)^{n+n'}` and that of :math:`K^{11}` and :math:`K^{22}` by
:math:`-(-1)^{n+n'}`, so one hemisphere gives half of the nonzero integrals.
Without mirror symmetry all four integrals are taken for every :math:`n,n'`
over the whole sphere, with :math:`a_{nn'}` in place of :math:`2a_{nn'}`.
All four vanish identically for :math:`m=0` except :math:`K^{12}` and
:math:`K^{21}`.  The regular integrals :math:`\mathrm{Rg}K^{pq}` are the same
with :math:`h_n` replaced by :math:`j_n(\xi_i)`.

With :math:`N=n_{\max}-n_0+1`, the :math:`2N\times2N` matrices

.. math::

   Q = \begin{pmatrix}
     i(\eta K^{21}-K^{12}) & -(\eta K^{11}+K^{22})\\
     -(\eta K^{22}+K^{11}) & i(K^{21}-\eta K^{12})
   \end{pmatrix}

and :math:`\mathrm{Rg}Q`, built likewise from :math:`\mathrm{Rg}K^{pq}`, give

.. math::

   T = -\,\mathrm{Rg}Q\;Q^{-1} = \begin{pmatrix} T^{11} & T^{12}\\ T^{21} & T^{22}\end{pmatrix},

with :math:`T^{pq}_{nn'}(m)` the blocks.  (The original multiplies :math:`Q` and
:math:`\mathrm{Rg}Q` by :math:`k^2`, which cancels in :math:`T`.)

.. _Sec T-matrix convergence:

Truncation and convergence
^^^^^^^^^^^^^^^^^^^^^^^^^^

The requested accuracy :math:`\delta` enters as :math:`\delta'=\delta/10`.
From the :math:`m=0` matrix,

.. math::

   q_{\mathrm{ext}} = \sum_{n=1}^{n_{\max}}(2n+1)\,\mathrm{Re}\bigl(T^{11}_{nn}+T^{22}_{nn}\bigr),\qquad
   q_{\mathrm{sca}} = \sum_{n=1}^{n_{\max}}(2n+1)\bigl(|T^{11}_{nn}|^2+|T^{22}_{nn}|^2\bigr).

Starting from
:math:`n_{\max}=\max(4,\lfloor x_v+4.05x_v^{1/3}\rfloor)`, :math:`x_v=kr_v`, and
:math:`N_G=N_D\,n_{\max}`, :math:`n_{\max}` increases by one until both
:math:`q` change relatively by at most :math:`\delta'`.  Then, at that
:math:`n_{\max}`, :math:`N_G` increases by one until the same holds.  The
convergence fails if :math:`\delta'` is below the double-precision epsilon,
where the changes are rounding; once an added node no longer reduces the change;
once the :math:`q` are not finite; and for random orientation once
:math:`\lfloor n_{\max}/2\rfloor+1`, a lower bound of :math:`n_1`, exceeds the
Wigner tables of the expansion.  :math:`N_D` is a parameter, 2 by default.

For random orientation, :math:`n_1` is the smallest
:math:`n\ge\lfloor n_{\max}/2\rfloor+1` whose terms in both :math:`q` are at
most :math:`\delta'` relative to the sums, or :math:`n_{\max}` if none is; it is
taken at the last :math:`n_{\max}` step.  Random orientation uses the orders
:math:`m\le n_1` and, in the scattering sums below, the degrees
:math:`n,n'\le n_1`.  Fixed orientation uses all orders and degrees up to
:math:`n_{\max}`.

.. _Sec T-matrix cross sections:

Orientation-averaged cross sections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   C_{\mathrm{ext}} &= -\frac{2\pi}{k^2}\sum_{m}(2-\delta_{m0})\sum_n
     \mathrm{Re}\bigl(T^{11}_{nn}(m)+T^{22}_{nn}(m)\bigr),\\
   C_{\mathrm{sca}} &= \frac{2\pi}{k^2}\sum_{m}(2-\delta_{m0})\sum_{n,n'}\sum_{p,q=1}^{2}
     \bigl|T^{pq}_{nn'}(m)\bigr|^2.

The factor :math:`2-\delta_{m0}` counts :math:`\pm m`.  For random orientation
the extinction sum runs over :math:`n\le n_{\max}` and the scattering sum over
:math:`n,n'\le n_1`, both for :math:`m\le n_1`.  Fixed orientation reports the
same averages with all orders and degrees.

.. _Sec T-matrix amplitude:

Amplitude matrix in a fixed orientation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The incident and scattered directions are :math:`(\theta_0,\phi_0)` and
:math:`(\theta,\phi)` in the laboratory frame, and the particle axis has the
Euler angles :math:`(\alpha,\beta)`.  In the particle frame a direction has
the components

.. math::

   k_x = \cos\beta\sin\theta\cos(\phi-\alpha)-\sin\beta\cos\theta,\quad
   k_y = \sin\theta\sin(\phi-\alpha),\quad
   k_z = \cos\theta\cos\beta+\sin\theta\sin\beta\cos(\phi-\alpha),

and the angles

.. math::

   \theta_p = \operatorname{atan2}\bigl(\sqrt{k_x^2+k_y^2},\;k_z\bigr),\qquad
   \phi_p = \operatorname{atan2}(k_y,\;k_x).

The angles are used as given, also along the particle axis, at the laboratory
poles and at :math:`\beta=\pi/2`.  :math:`\theta_p` from the components is
accurate to rounding next to the axis, where :math:`\arccos k_z` would turn the
rounding of :math:`k_z` into an error of order :math:`\sqrt{\epsilon}`, and
neither angle divides by :math:`\sin\theta_p`.  On the axis, :math:`\phi_p` is
arbitrary and :math:`S` does not depend on it, since it enters both :math:`A_p`
and :math:`\Delta\phi`.  At a laboratory pole, the basis is the limit along the
given azimuth.

With :math:`\Delta\phi=\phi_p-\phi_{0,p}`,
:math:`c_{nn'}=i^{\,n'-n-1}\sqrt{(2n+1)(2n'+1)/(n(n+1)n'(n'+1))}`, the
functions :math:`\pi_n,\tau_n` of order :math:`m` at :math:`\theta_p` and
:math:`\pi^0_{n'},\tau^0_{n'}` at :math:`\theta_{0,p}`, and :math:`T^{pq}` the
blocks of order :math:`m` at :math:`(n,n')`, the particle-frame amplitude
matrix is

.. math::

   kS^{p}_{vv} &= \sum_{n,n'} c_{nn'}\tau_n\tau^0_{n'}T^{22}_{nn'}(0)
     + \sum_{m\ge1}2\cos(m\Delta\phi)\sum_{n,n'}c_{nn'}
       \bigl[T^{11}\pi_n\pi^0_{n'}+T^{21}\tau_n\pi^0_{n'}+T^{12}\pi_n\tau^0_{n'}+T^{22}\tau_n\tau^0_{n'}\bigr],\\
   kS^{p}_{vh} &= \sum_{m\ge1}2\sin(m\Delta\phi)\sum_{n,n'}c_{nn'}
       \bigl[T^{11}\pi_n\tau^0_{n'}+T^{21}\tau_n\tau^0_{n'}+T^{12}\pi_n\pi^0_{n'}+T^{22}\tau_n\pi^0_{n'}\bigr],\\
   kS^{p}_{hv} &= -\sum_{m\ge1}2\sin(m\Delta\phi)\sum_{n,n'}c_{nn'}
       \bigl[T^{11}\tau_n\pi^0_{n'}+T^{21}\pi_n\pi^0_{n'}+T^{12}\tau_n\tau^0_{n'}+T^{22}\pi_n\tau^0_{n'}\bigr],\\
   kS^{p}_{hh} &= \sum_{n,n'} c_{nn'}\tau_n\tau^0_{n'}T^{11}_{nn'}(0)
     + \sum_{m\ge1}2\cos(m\Delta\phi)\sum_{n,n'}c_{nn'}
       \bigl[T^{11}\tau_n\tau^0_{n'}+T^{21}\pi_n\tau^0_{n'}+T^{12}\tau_n\pi^0_{n'}+T^{22}\pi_n\pi^0_{n'}\bigr],

with :math:`n,n'\ge\max(m,1)`.  In matrix form, with the block
:math:`T_{nn'}(m)=\begin{pmatrix}T^{11}&T^{12}\\T^{21}&T^{22}\end{pmatrix}`
and :math:`P_n=\begin{pmatrix}\pi_n&\tau_n\\\tau_n&\pi_n\end{pmatrix}`, and
:math:`P^0_{n'}` likewise of the incident direction, order :math:`m` contributes
:math:`M_m=\sum_{n,n'}c_{nn'}P_nT_{nn'}(m)P^0_{n'}`: its co-polar elements times
:math:`(2-\delta_{m0})\cos(m\Delta\phi)`, and its cross-polar elements times
:math:`(2-\delta_{m0})\sin(m\Delta\phi)` above and its negative below the
diagonal.  The laboratory amplitude matrix is
:math:`S=R_s^{-1}S^pR_0`, where :math:`R=A_p\,B\,A_l` for the incident
(:math:`R_0`) and scattered (:math:`R_s`) direction, with

.. math::

   B = \begin{pmatrix}
     \cos\alpha\cos\beta & \sin\alpha\cos\beta & -\sin\beta\\
     -\sin\alpha & \cos\alpha & 0\\
     \cos\alpha\sin\beta & \sin\alpha\sin\beta & \cos\beta
   \end{pmatrix},\quad
   A_l = \begin{pmatrix}
     \cos\theta\cos\phi & -\sin\phi\\
     \cos\theta\sin\phi & \cos\phi\\
     -\sin\theta & 0
   \end{pmatrix},\quad
   A_p = \begin{pmatrix}
     \cos\theta_p\cos\phi_p & \cos\theta_p\sin\phi_p & -\sin\theta_p\\
     -\sin\phi_p & \cos\phi_p & 0
   \end{pmatrix}.

:math:`S` has the units of length.  ARTS's meridional basis of a ray is
:math:`(\hat v,\hat h)=(-\hat\theta,-\hat\phi)` in this frame, so :math:`S`
is the Jones matrix between ARTS's bases of the two rays, and the phase matrix
is its Mueller matrix (:doc:`concept.polarization`).

.. _Sec T-matrix aro:

Azimuthally random orientation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ARTS's propagation direction of zenith angle :math:`\vartheta` and azimuth
:math:`\varphi`, clockwise from north in the east-north-up frame, has the basis
vector :math:`\hat h=(\cos\varphi,-\sin\varphi,0)`
(:doc:`concept.polarization`).  It is the laboratory direction above with
:math:`\theta=\vartheta` and :math:`\phi=90^\circ-\varphi`, and the bases
coincide.  The data of an azimuthally random particle depend on the azimuths
only through the scattered minus the incident one, :math:`\Delta\varphi`, so
the amplitude matrix is evaluated at :math:`\phi_0=90^\circ` and
:math:`\phi=90^\circ-\Delta\varphi`.

The particle's axis has the tilts :math:`\beta_j` from the vertical with the
weights :math:`w_j`, normalized to one, and every azimuth :math:`\alpha` alike.
A quantity :math:`X` of the particle, the phase matrix
:math:`Z=M(S)` or the extinction matrix and absorption vector below, is
averaged as

.. math::

   \langle X\rangle = \sum_j w_j\,\frac{1}{N_\alpha}\sum_{k=0}^{N_\alpha-1}X(\alpha_k,\beta_j),\qquad
   \alpha_k = \frac{2\pi k}{N_\alpha},\quad N_\alpha=4n_{\max}+1.

Rotating the particle by :math:`\alpha` about the laboratory :math:`z` axis
multiplies the laboratory T-matrix element of the orders :math:`m` and
:math:`m'` by :math:`e^{-i(m-m')\alpha}`, with :math:`|m|,|m'|\le n_{\max}`.
:math:`S` is therefore a trigonometric polynomial of degree :math:`2n_{\max}`
in :math:`\alpha`, and :math:`Z` of degree :math:`4n_{\max}`.  The trapezoidal
rule of :math:`N` nodes integrates trigonometric polynomials of degree below
:math:`N` exactly, so the average is exact for the truncated T-matrix.  A
vertical axis, :math:`\beta=0` or :math:`180^\circ`, is invariant to
:math:`\alpha`, and one node is exact.

The ensemble is mirror symmetric in every vertical plane: the mirror image of
a body of revolution is the same body with its axis mirrored, and that axis
has the same tilt.  The averaged forward amplitude
:math:`S(\hat k\leftarrow\hat k)` is therefore diagonal, and the extinction
matrix has the elements

.. math::

   K_{11} = \lambda\,\mathrm{Im}(S_{vv}+S_{hh}),\qquad
   K_{12} = \lambda\,\mathrm{Im}(S_{vv}-S_{hh}),\qquad
   K_{34} = \lambda\,\mathrm{Re}(S_{vv}-S_{hh}),

with :math:`K_{21}=K_{12}` and :math:`K_{43}=-K_{34}`.  These are the linear
term of the coherent transmission,
:math:`K=\tfrac{\lambda}{2}\bigl[M(I-iS)-M(I+iS)\bigr]` with :math:`M` the
Mueller matrix in ARTS's Stokes convention, which fixes the sign of
:math:`K_{34}`.

The absorption vector is the emission that keeps unpolarized isotropic
radiation in equilibrium: it balances the extinction less the scattering into
the direction :math:`\hat k` from all others,

.. math::

   \boldsymbol a(\hat k) = K(\hat k)\,\boldsymbol e_1-\int_{4\pi}Z(\hat k\leftarrow\hat k')\,\boldsymbol e_1\,d\Omega',
   \qquad \boldsymbol e_1=(1,0,0,0)^T,

of which only :math:`a_I` and :math:`a_Q` are nonzero, by the mirror symmetry.  As a
function of :math:`\hat k'`, :math:`Z\boldsymbol e_1` sums over the incident
polarizations and has spherical-harmonic degree at most :math:`2n_{\max}`:
the :math:`(n_{\max}+1)`-point Gauss-Legendre rule in :math:`\cos\theta'` and
the :math:`(2n_{\max}+1)`-point trapezoidal rule in :math:`\phi'` integrate it
exactly.  For a sphere, :math:`a_I` is the absorption cross section.

.. _Sec T-matrix expansion:

Expansion coefficients in random orientation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The scattering matrix of the randomly oriented particle is expanded in
generalized spherical functions (:ref:`Sec T-matrix matrix`).  Its coefficients
follow from the T-matrix analytically.  :math:`C(j_1m_1j_2m_2|JM)` is the
Clebsch-Gordan coefficient with the Condon-Shortley phase.  All degrees and
orders below are at most :math:`n_1` unless stated.

For :math:`0\le m\le\min(n,n')` and :math:`k=1,2`, with the blocks of order
:math:`m`,

.. math::

   T^{(1)}_{nn'}(\pm m) &= T^{11}_{nn'}\pm T^{12}_{nn'}\pm T^{21}_{nn'}+T^{22}_{nn'},\\
   T^{(2)}_{nn'}(\pm m) &= T^{11}_{nn'}\pm T^{12}_{nn'}\mp T^{21}_{nn'}-T^{22}_{nn'},

the first column of :math:`UT_{nn'}U` for :math:`+m` and of
:math:`UZT_{nn'}ZU` for :math:`-m`, with :math:`U=\begin{pmatrix}1&1\\1&-1\end{pmatrix}`
and :math:`Z=\operatorname{diag}(1,-1)`.

For each :math:`n`, each :math:`n_s=0,\ldots,n_1+n` and
:math:`n'=\max(1,|n-n_s|),\ldots,\min(n_1,n+n_s)`,

.. math::

   A^{(k)}_{nn_sn'} = \frac{i^{\,n'-n}}{\sqrt{2n'+1}}
     \sum_{m=0}^{\min(n,n')}C(n\,m\,n_s\,0|n'\,m)
     \Bigl[T^{(k)}_{nn'}(m)+(1-\delta_{m0})(-1)^{n+n_s+n'}T^{(k)}_{nn'}(-m)\Bigr],

and for :math:`\max(1-n_s,-n)\le m\le\min(n_s+1,n)`, over the same :math:`n'`,

.. math::

   B^{(k)}_{n_smn} = \sum_{n'} C(n\,m\,n_s\,1{-}m|n'\,1)\,A^{(k)}_{nn_sn'}.

For each pair :math:`n,\nu`, with sums over
:math:`n_s=|m-1|,\ldots,n_1+\min(n,\nu)`,

.. math::

   D^{(k)}_{m\nu n} &= \sum_{n_s}(2n_s+1)\,\mathrm{Re}\bigl(B^{(k)}_{n_smn}B^{(k)*}_{n_sm\nu}\bigr),
     &&|m|\le\min(n,\nu),\\
   D^{(2+k)}_{m\nu n} &= \sum_{n_s}(2n_s+1)\,\mathrm{Re}\bigl(B^{(k)}_{n_smn}B^{(k)*}_{n_s,2-m,\nu}\bigr),
     &&\max(-n,2-\nu)\le m\le\min(n,\nu+2),\\
   D^{(5)}_{m\nu n} &= \sum_{n_s}(2n_s+1)\,B^{(2)}_{n_smn}B^{(1)*}_{n_s,2-m,\nu},
     &&\text{same } m.

Then for :math:`l=0,1,\ldots,2n_1`, summing over :math:`n` and
:math:`\nu=\max(1,|n-l|),\ldots,\min(n_1,n+l)`, with
:math:`\sigma=(-1)^{n+l+\nu}`, :math:`f=\sqrt{(2n+1)/(2\nu+1)}`,
:math:`c_1=f\,C(n\,1\,l\,0|\nu\,1)` and :math:`c_2=f\,C(n\,{-1}\,l\,2|\nu\,1)`,

.. math::

   g_1 &= \sum c_1\sum_{m}C(n\,m\,l\,0|\nu\,m)\,D^{(1)}_{m\nu n}, &
   g_2 &= \sum c_1\,\sigma\sum_{m}C(n\,m\,l\,0|\nu\,m)\,D^{(2)}_{m\nu n},\\
   g_3 &= \sum c_2\sum_{m}C(n\,{-m}\,l\,2|\nu\,2{-}m)\,D^{(3)}_{m\nu n}, &
   g_4 &= \sum c_2\,\sigma\sum_{m}C(n\,{-m}\,l\,2|\nu\,2{-}m)\,D^{(4)}_{m\nu n},\\
   g_5 &= -\sum c_1\sum_{m}C(n\,{-m}\,l\,2|\nu\,2{-}m)\,D^{(5)}_{m\nu n},

over the ranges of :math:`m` above; :math:`g_3,g_4,g_5` are zero for
:math:`l<2`.  The coefficients of :math:`m<0` in :math:`g_1` and :math:`g_2`
follow from :math:`C(n\,{-m}\,l\,0|\nu\,{-m})=\sigma\,C(n\,m\,l\,0|\nu\,m)`,
as in the original, and :math:`c_1` and :math:`c_2` are the :math:`m=1` terms
of their sums.  With :math:`c_l=(2l+1)\,\lambda^2/(4\pi C_{\mathrm{sca}})`,

.. math::

   \alpha_1^l = c_l(g_1+g_2),\quad \alpha_4^l = c_l(g_1-g_2),\quad
   \alpha_2^l = c_l(g_3+g_4),\quad \alpha_3^l = c_l(g_3-g_4),\quad
   \beta_1^l = 2c_l\,\mathrm{Re}\,g_5,\quad \beta_2^l = 2c_l\,\mathrm{Im}\,g_5.

The series has all :math:`2n_1+1` coefficients, :math:`l=0,\ldots,2n_1`, the
exact series of the truncated T-matrix.  Normalization gives :math:`\alpha_1^0=1`, and the asymmetry parameter is
:math:`\alpha_1^1/3`.

.. _Sec T-matrix matrix:

Scattering matrix from the expansion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the scattering angles :math:`\Theta_j=\pi(j-1)/(N_\Theta-1)`, with
:math:`u=\cos\Theta`,

.. math::

   F_{11} = \sum_l\alpha_1^lP^l_{00}(u),\quad F_{44} = \sum_l\alpha_4^lP^l_{00}(u),\quad
   F_{12} = \sum_{l\ge2}\beta_1^lP^l_{02}(u),\quad F_{34} = \sum_{l\ge2}\beta_2^lP^l_{02}(u),

.. math::

   F_{22}\pm F_{33} = \sum_{l\ge2}(\alpha_2^l\pm\alpha_3^l)\,P^l_{2,\pm2}(u),

with :math:`P^l_{00}=P_l` and the generalized spherical functions from
:math:`P^2_{22}=\tfrac14(1+u)^2`, :math:`P^2_{2,-2}=\tfrac14(1-u)^2`,
:math:`P^2_{02}=\tfrac{\sqrt6}{4}(u^2-1)` by

.. math::

   P^{l+1}_{2,\pm2} &= \frac{(2l+1)\bigl(l(l+1)u\mp4\bigr)P^l_{2,\pm2}-(l+1)(l^2-4)P^{l-1}_{2,\pm2}}{l\bigl((l+1)^2-4\bigr)},\\
   P^{l+1}_{02} &= \frac{(2l+1)\,u\,P^l_{02}-\sqrt{l^2-4}\,P^{l-1}_{02}}{\sqrt{(l+1)^2-4}},

starting from zero at :math:`l=1`.  These are Mishchenko's elements; in ARTS's
Stokes convention :math:`F_{34}=-\sum_{l\ge2}\beta_2^lP^l_{02}(u)`.

.. _Sec T-matrix size distribution:

Size distributions
^^^^^^^^^^^^^^^^^^

A size distribution :math:`n(r)` over :math:`[r_1,r_2]` is integrated with the
Gauss-Legendre rule of :math:`N_k` points, two of them for the modified power
law (on :math:`[0,r_1]` and :math:`[r_1,r_2]`).  The weights are multiplied by
:math:`n(r)` and normalized to one.  With the parameters :math:`a`, :math:`b`
and :math:`\gamma`, the distributions are

.. math::

   \text{modified gamma:}&\quad n(r)\propto r^{a}\exp\bigl(-(a/\gamma)(r/b)^{\gamma}\bigr),\\
   \text{log-normal:}&\quad n(r)\propto r^{-1}\exp\bigl(-\ln^2(r/a)/(2b)\bigr),\\
   \text{power law:}&\quad n(r)\propto r^{-3},\\
   \text{gamma:}&\quad n(r)\propto r^{(1-3b)/b}\exp\bigl(-r/(ab)\bigr),\\
   \text{modified power law:}&\quad n(r)\propto 1\ (r\le r_1),\quad (r/r_1)^{b}\ (r>r_1).

For the power law of :cite:t:`Hansen_Travis:74`, :math:`a` and :math:`b` are the
effective radius and variance, and :math:`r_1` solves
:math:`(r_2-r_1)/\ln(r_2/r_1)=a` with :math:`r_2=2a(1+b)-r_1` on
:math:`(0,a]`.  The left side is :math:`-a` in the limit :math:`r_1\to0` and
positive at :math:`r_1=a` for every :math:`b>0`, so the root exists; it is about
:math:`2a(1+b)e^{-2(1+b)}`.  Boost's TOMS 748 finds it.  It at least halves its
bracket in each pass of at most four evaluations, so any root above the
subnormal range is resolved within :math:`2+4(e_{\max}-e_{\min}+p)` evaluations
of the double format's exponent range and precision.  The distribution fails if
:math:`1+2b` rounds to 1, or if the root is not resolved within them.  For the
other distributions :math:`r_1` and
:math:`r_2` are given multiples of :math:`a`.  The effective radius and variance of
the distribution are

.. math::

   r_{\mathrm{eff}} = \frac{\sum_i r_i^3w_i}{\sum_i r_i^2w_i},\qquad
   v_{\mathrm{eff}} = \frac{\sum_i(r_i-r_{\mathrm{eff}})^2r_i^2w_i}{r_{\mathrm{eff}}^2\sum_i r_i^2w_i}.

The node radii :math:`r_i` are of the kind given (volume- or
surface-area-equivalent); the particle at node :math:`i` has the
volume-equivalent radius :math:`q\,r_i` (:ref:`Sec T-matrix shapes`).  The
averages are

.. math::

   C_{\mathrm{sca}} = \sum_iw_iC_{\mathrm{sca},i},\quad
   C_{\mathrm{ext}} = \sum_iw_iC_{\mathrm{ext},i},\quad
   \alpha^l = \frac{1}{C_{\mathrm{sca}}}\sum_iw_iC_{\mathrm{sca},i}\,\alpha^l_i,

for every coefficient, the single-scattering albedo is
:math:`C_{\mathrm{sca}}/C_{\mathrm{ext}}`, and the averaged coefficients give the
scattering matrix as above.
