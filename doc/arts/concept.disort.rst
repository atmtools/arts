.. _Sec DISORT:

DISORT and VDISORT core solvers
################################

ARTS contains two closely related plane-parallel discrete-ordinate solvers in
``src/core/disort-cpp``:

* ``disort.cpp`` solves the scalar radiative-transfer equation; and
* ``vdisort.cpp`` solves its polarized, Stokes-vector form.

This page documents the mathematics and numerical conventions of those core
solvers.  The workspace methods which prepare their inputs and expose their
outputs are documented separately.

Both implementations follow the traditional DISORT separation into azimuthal
Fourier modes, a layerwise eigenproblem, and one banded boundary-value solve per
mode described by `Stamnes et al. (1988)
<https://doi.org/10.1364/AO.27.002502>`_.  Their implementation and validation
draw on the original DISORT 4.0.99 code, cDISORT, and `PythonicDISORT
<https://doi.org/10.21105/joss.06442>`_.  The vector mathematics follows the
VSIDORT full-eigenproblem formulation of `Lin et al. (2022)
<https://doi.org/10.3389/frsen.2022.880768>`_.  The remainder of this page
states the equations and conventions used by the ARTS solvers directly rather
than describing them through those source packages.

Geometry and sign conventions
*****************************

The atmosphere consists of homogeneous plane-parallel layers.  Optical depth
:math:`\tau` is zero at the top and increases downwards.  Direction cosines are
defined by

.. math::

   \mu = \cos\theta,

with :math:`\mu>0` directed upwards and :math:`\mu<0` directed downwards.  A
direct beam with positive :math:`\mu_0` therefore travels in direction
:math:`-\mu_0`.  Azimuth :math:`\phi` is in radians and is reconstructed using
:math:`\phi_0-\phi`, where :math:`\phi_0` is the beam azimuth.

For an even number :math:`N_q` of streams, the solver uses :math:`N=N_q/2`
positive Gauss--Legendre nodes and their negatives.  Internally the positive
streams precede the corresponding negative streams.  The positive-node
weights :math:`w_j` are reused for both hemispheres.

More explicitly, the direction and weight arrays are

.. math::

   \boldsymbol\mu
     &= (\mu_1,\ldots,\mu_N,-\mu_1,\ldots,-\mu_N),
        \qquad \mu_i>0,\\
   \overline{\boldsymbol w}
     &= (w_1,\ldots,w_N,w_1,\ldots,w_N).

Mathematical indices below start at one.  The corresponding C++ stream indices
are :math:`i=0,\ldots,N-1` for :math:`+\mu_i` and :math:`N+i` for
:math:`-\mu_i`.  With :math:`L` layers, define the physical layer edges by

.. math::

   0=t_0<t_1<\cdots<t_L,

where ``tau_arr[l]`` stores :math:`t_{l+1}`.  Scalar DISORT uses the
delta-M-scaled coordinate :math:`t^*` inside its transport eigensystem;
VDISORT uses the optical-depth coordinate supplied to it.

The geometry is deliberately plane-parallel.  Neither core implements the
pseudo-spherical direct-beam approximation.  Curved ray
paths belong in the ARTS geometry layer rather than in a solver tied to a
single planetary radius.

Scalar radiative transfer
*************************

For scalar radiance :math:`I(\tau,\mu,\phi)`, the equation solved in a layer is

.. math::

   \mu\frac{\partial I}{\partial\tau}
   = I
   - \frac{\omega}{4\pi}\int_{4\pi}
       P(\Omega,\Omega') I(\Omega')\,\mathrm{d}\Omega'
   - S_{\mathrm{dir}} - S_{\mathrm{int}},

where :math:`\omega` is the single-scattering albedo,
:math:`P` is normalized to an integral of :math:`4\pi`, and the last two terms
are the direct-beam and internal sources.  The phase function is represented by
Legendre moments

.. math::

   P(\cos\Theta)
   = \sum_{n=0}^{N_{\mathrm{leg}}-1}(2n+1)\chi_n P_n(\cos\Theta),
   \qquad \chi_0=1.

The azimuthal dependence is expanded as

.. math::

   I(\tau,\mu,\phi)
   = \sum_{m=0}^{M-1} I^m(\tau,\mu)
       \cos\!\left[m(\phi_0-\phi)\right].

After applying Gauss--Legendre quadrature, each Fourier mode has the form

.. math::

   \frac{\mathrm{d}\boldsymbol I^m}{\mathrm{d}\tau}
   = \boldsymbol A^m\boldsymbol I^m
     + \boldsymbol q^m(\tau).

Here :math:`\boldsymbol q^m` is the forcing *after* division by the signed
stream cosine.  Its signs therefore depend on the hemisphere.  The physical
direct-beam and internal sources occur with a minus sign in the transfer
equation written above; diffuse scattering has already been included in
:math:`\boldsymbol A^m`.

Scalar stream matrix
====================

Let

.. math::

   \boldsymbol I^m=
   \begin{bmatrix}\boldsymbol I^{m,+}\\
                   \boldsymbol I^{m,-}\end{bmatrix},\qquad
   \boldsymbol M=\operatorname{diag}(\mu_1,\ldots,\mu_N),\qquad
   \boldsymbol W=\operatorname{diag}(w_1,\ldots,w_N).

For layer :math:`\ell`, define the same- and opposite-hemisphere scattering
blocks

.. math::

   D^{m,+}_{ij}
   &= \frac{\omega_\ell^*}{2}
      \sum_{n=m}^{N_{\mathrm{leg}}-1}
      (2n+1)\chi^*_{\ell,n}\frac{(n-m)!}{(n+m)!}
      P_n^m(\mu_i)P_n^m(\mu_j),\\
   D^{m,-}_{ij}
   &= \frac{\omega_\ell^*}{2}
      \sum_{n=m}^{N_{\mathrm{leg}}-1}
      (2n+1)\chi^*_{\ell,n}\frac{(n-m)!}{(n+m)!}
      P_n^m(\mu_i)P_n^m(-\mu_j).

Stars denote delta-M-scaled quantities.  These matrices are ``D_pos`` and
``D_neg``.  The factorial ratio is explicit because the implementation uses
ordinary rather than semi-normalized associated Legendre functions.

The conceptual left and right matrices of the full generalized eigenproblem
are

.. math::

   \underbrace{\left[
   \boldsymbol I_{2N}-
   \begin{pmatrix}
     \boldsymbol D^{m,+}\boldsymbol W &
     \boldsymbol D^{m,-}\boldsymbol W\\
     \boldsymbol D^{m,-}\boldsymbol W &
     \boldsymbol D^{m,+}\boldsymbol W
   \end{pmatrix}\right]}_{\boldsymbol L_{\mathrm{eig}}^m}
   \boldsymbol g
   = k\,
   \underbrace{\begin{pmatrix}
     \boldsymbol M&0\\0&-\boldsymbol M
   \end{pmatrix}}_{\boldsymbol R_{\mathrm{eig}}}
   \boldsymbol g.

CPP-DISORT does not pass this generalized pair to LAPACK.  It first defines

.. math::

   \boldsymbol\alpha
     &=\boldsymbol M^{-1}
       (\boldsymbol D^{m,+}\boldsymbol W-\boldsymbol I_N),\\
   \boldsymbol\beta
     &=\boldsymbol M^{-1}\boldsymbol D^{m,-}\boldsymbol W,

so that premultiplication by :math:`\boldsymbol R_{\mathrm{eig}}^{-1}`
gives the ordinary :math:`2N\times2N` matrix

.. math::

   \boldsymbol A^m=
   \begin{pmatrix}
     -\boldsymbol\alpha&-\boldsymbol\beta\\
      \boldsymbol\beta& \boldsymbol\alpha
   \end{pmatrix},\qquad
   \boldsymbol A^m\boldsymbol g=k\boldsymbol g.

Scalar reduced eigenproblem
===========================

The full matrix is also not explicitly diagonalized.  If

.. math::

   \boldsymbol x=\boldsymbol g^++\boldsymbol g^-,\qquad
   \boldsymbol y=\boldsymbol g^--\boldsymbol g^+,

then

.. math::

   (\boldsymbol\alpha+\boldsymbol\beta)\boldsymbol x
      &=k\boldsymbol y,\\
   (\boldsymbol\alpha-\boldsymbol\beta)\boldsymbol y
      &=k\boldsymbol x.

Consequently the actual eigensystem passed to the real eigensolver is only
:math:`N\times N`:

.. math::

   (\boldsymbol\alpha-\boldsymbol\beta)
   (\boldsymbol\alpha+\boldsymbol\beta)\boldsymbol X
   =\boldsymbol X\operatorname{diag}(\boldsymbol\kappa^2).

Reversing the product and using the complementary sum/difference eigenvector
is an equivalent reduction and recovers the same full spectrum.  The order
above is the one used by this implementation.

CPP-DISORT stores the paired propagation constants as

.. math::

   \boldsymbol K=
   (-\sqrt{|\kappa_1^2|},\ldots,-\sqrt{|\kappa_N^2|},
     +\sqrt{|\kappa_1^2|},\ldots,+\sqrt{|\kappa_N^2|}).

The absolute value is the literal implementation.  For a physically valid
scalar system the squared eigenvalues are expected to be non-negative; taking
the absolute value prevents a small negative round-off value from producing a
NaN, but can also conceal a genuinely invalid reduced eigensystem.

With

.. math::

   \boldsymbol Z=(\boldsymbol\alpha+\boldsymbol\beta)
      \boldsymbol X\operatorname{diag}(\boldsymbol\kappa^{-1}),

the full eigenvector matrix, with negative propagation constants first, is

.. math::

   \boldsymbol G=\frac12
   \begin{pmatrix}
     \boldsymbol X+\boldsymbol Z&\boldsymbol X-\boldsymbol Z\\
     \boldsymbol X-\boldsymbol Z&\boldsymbol X+\boldsymbol Z
   \end{pmatrix}.

In the implementation, ``apb`` and ``amb`` represent
:math:`\boldsymbol\alpha+\boldsymbol\beta` and
:math:`\boldsymbol\alpha-\boldsymbol\beta`, respectively.

Scalar layer solution
=====================

The homogeneous part in layer :math:`\ell` is

.. math::

   \boldsymbol I^m_{\ell,\mathrm{hom}}(\tau^*)
   &=\sum_{e=1}^{N}C_{\ell,e}\boldsymbol G_{\ell,e}
       \exp[K_{\ell,e}(\tau^*-t_\ell^*)]\\
   &\quad+
     \sum_{e=N+1}^{2N}C_{\ell,e}\boldsymbol G_{\ell,e}
       \exp[K_{\ell,e}(\tau^*-t_{\ell+1}^*)].

Thus negative propagation constants are anchored at the layer top and
positive constants at the layer bottom.  The direct-beam particular solution
has the form

.. math::

   \boldsymbol I^m_{\ell,\mathrm b}(\tau^*)
      =\boldsymbol B_\ell^m\exp(-\tau^*/\mu_0).

If :math:`\widetilde{\boldsymbol X}^m` denotes its signed-stream forcing
amplitude in the standard ODE, then

.. math::

   (\boldsymbol A^m+\mu_0^{-1}\boldsymbol I_{2N})
      \boldsymbol B_\ell^m=-\widetilde{\boldsymbol X}^m.

CPP-DISORT evaluates this solve in the eigenbasis.  The isotropic polynomial
particular solution exists only for :math:`m=0`; it is likewise transformed
through :math:`\boldsymbol G^{-1}` and evaluated by a backwards polynomial
recurrence.  The complete layer field is the sum of homogeneous, beam, and
polynomial parts.

Polarized radiative transfer
****************************

VDISORT transports the ARTS Stokes vector

.. math::

   \boldsymbol I = [I,Q,U,V]^{\mathsf T},

not the alternative :math:`[I_v,I_h,U,V]^{\mathsf T}` representation.  With
the linear-polarization convention used here, the two are related by

.. math::

   \begin{bmatrix}I\\Q\end{bmatrix}
      =\begin{pmatrix}1&1\\1&-1\end{pmatrix}
       \begin{bmatrix}I_v\\I_h\end{bmatrix},\qquad
   \begin{bmatrix}I_v\\I_h\end{bmatrix}
      =\frac12\begin{pmatrix}1&1\\1&-1\end{pmatrix}
       \begin{bmatrix}I\\Q\end{bmatrix}.

Mueller matrices from a derivation using :math:`I_v,I_h` must be transformed
on both sides by this change of basis before they are identified with the ARTS
matrices.  The vector equation in the ARTS basis is

.. math::

   \mu\frac{\partial\boldsymbol I}{\partial\tau}
   = \boldsymbol I
   - \frac{\omega}{4\pi}\int_{4\pi}
       \mathbf P(\Omega,\Omega')\boldsymbol I(\Omega')
       \,\mathrm{d}\Omega'
   - \boldsymbol S_{\mathrm{dir}}
   - \boldsymbol S_{\mathrm{int}},

where :math:`\mathbf P` is a Mueller matrix expressed in a consistent Stokes
reference frame.

The Fourier expansion of a Mueller problem couples ordinary cosine and sine
coefficients.  VDISORT therefore stores two *combined* systems, conventionally
indexed by :math:`\alpha=0` and :math:`\alpha=1`, for every Fourier order.  At
:math:`m=0`, the first carries the :math:`I,Q` block and the second the
:math:`U,V` block.  At :math:`m>0`, the combined matrices mix the ordinary
cosine and sine Mueller coefficients according to the block layout stated
below.  The physical field is reconstructed schematically as

.. math::

   \begin{aligned}
   \begin{pmatrix} I \\ Q \end{pmatrix}
     &= \sum_m\left(\boldsymbol C^m_{IQ}\cos m\Delta\phi
                    +\boldsymbol S^m_{IQ}\sin m\Delta\phi\right),\\
   \begin{pmatrix} U \\ V \end{pmatrix}
     &= \sum_m\left(\boldsymbol S^m_{UV}\cos m\Delta\phi
                    +\boldsymbol C^m_{UV}\sin m\Delta\phi\right),
   \end{aligned}

with :math:`\Delta\phi=\phi_0-\phi`.  The swapped roles in the second line are
important; treating every Stokes component like scalar intensity gives the
wrong signs and azimuthal behavior.

The exact combined Mueller layout can be seen by partitioning Stokes space
into :math:`A=(I,Q)` and :math:`B=(U,V)`.  If
:math:`\boldsymbol P_{\cos}` and :math:`\boldsymbol P_{\sin}` are written in
:math:`2\times2` Stokes blocks, then for :math:`m>0` the two matrices supplied
to the independent transport systems are

.. math::

   \boldsymbol P^{c,m}
   &=\begin{pmatrix}
      (P_{\cos})_{AA}&-(P_{\sin})_{AB}\\
      (P_{\sin})_{BA}& (P_{\cos})_{BB}
     \end{pmatrix},\\
   \boldsymbol P^{s,m}
   &=\begin{pmatrix}
      (P_{\cos})_{AA}& (P_{\sin})_{AB}\\
     -(P_{\sin})_{BA}& (P_{\cos})_{BB}
     \end{pmatrix}.

For :math:`m=0`, they reduce to

.. math::

   \boldsymbol P^{c,0}
      =\begin{pmatrix}(P_{\cos})_{AA}&0\\0&0\end{pmatrix},\qquad
   \boldsymbol P^{s,0}
      =\begin{pmatrix}0&0\\0&(P_{\cos})_{BB}\end{pmatrix}.

Thus :math:`c` and :math:`s` label Lin's *combined systems*.  They must not be
interpreted as four independent ordinary cosine and sine Stokes vectors.

Scattering data produced in a complex spectral representation are split as

.. math::

   \mathbf P_{\cos}=\Re\mathbf P_{\mathbb C},\qquad
   \mathbf P_{\sin}=-\Im\mathbf P_{\mathbb C}.

No factor of two is introduced by this split.  Fourier normalization is handled
when the modes are combined and synthesized.  This interpretation applies to
the phase-matrix output, rather than changing the scattering-species data
model.

VDISORT stream and Stokes matrix
================================

Let :math:`D=4N_q=8N`.  Stream :math:`i` and Stokes component
:math:`s\in\{I,Q,U,V\}` are flattened in stream-major order,

.. math::

   \rho(i,s)=4i+s.

For example, all four components at :math:`+\mu_1` precede all four at
:math:`+\mu_2`; the downward streams start at flattened index :math:`4N`.
Define

.. math::

   \widehat{\boldsymbol M}
      &=\operatorname{diag}(\mu_1,\ldots,\mu_N,
                            -\mu_1,\ldots,-\mu_N)\otimes\boldsymbol I_4,\\
   \widehat{\boldsymbol W}
      &=\operatorname{diag}(w_1,\ldots,w_N,
                             w_1,\ldots,w_N)\otimes\boldsymbol I_4.

For a fixed combined system :math:`a\in\{c,s\}`, Fourier order :math:`m`,
and layer :math:`\ell`, let :math:`\boldsymbol{\mathcal P}^{am\ell}` be the
:math:`D\times D` block matrix whose :math:`(i,j)` block is the Mueller matrix
:math:`\boldsymbol P^{am\ell}_{ij}`.  The conceptual generalized eigenproblem is

.. math::

   \underbrace{\left[
      \boldsymbol I_D-
      \frac{\omega_\ell}{2}
      \boldsymbol{\mathcal P}^{am\ell}\widehat{\boldsymbol W}
   \right]}_{\boldsymbol L_{\mathrm{eig}}^{am\ell}}
   \boldsymbol g
   =k\,
   \underbrace{\widehat{\boldsymbol M}}_{\boldsymbol R_{\mathrm{eig}}}
   \boldsymbol g.

As in the scalar core, this pair is not sent to a generalized eigensolver.
VDISORT explicitly forms

.. math::

   \boldsymbol A^{am\ell}
      =\widehat{\boldsymbol M}^{-1}
       \boldsymbol L_{\mathrm{eig}}^{am\ell},

whose individual entries are

.. math::

   A^{am\ell}_{\rho(i,s_o),\rho(j,s_i)}
   =\frac{\delta_{ij}\delta_{s_os_i}}{\mu_i}
    -\frac{\omega_\ell}{2\mu_i}\,
       \overline w_j
       \left[P^{am\ell}_{ij}\right]_{s_os_i}.

Here :math:`\mu_i` is signed and :math:`\overline w_j` is the repeated weight
defined with the stream array above (``W[j % N]`` in C++).  A Mueller block
maps incident to outgoing Stokes components, so :math:`s_o` is the matrix row
and :math:`s_i` its column.

VDISORT full eigenproblem
=========================

The implementation solves the ordinary, generally non-symmetric problem

.. math::

   \boldsymbol A^{am\ell}\boldsymbol G^{am\ell}
      =\boldsymbol G^{am\ell}\operatorname{diag}(\boldsymbol K^{am\ell})

as a full :math:`D\times D` complex eigensystem.  There are two such
eigensystems for each :math:`(m,r)`, one for each combined system; there is not
one :math:`2D\times2D` system.  The input matrix is real, but general Mueller
coupling makes it non-symmetric, so complex conjugate eigenpairs are retained.
This differs from an explicitly real conjugate-pair basis, but reconstructs the
same real physical field when those pairs cancel.

Eigenvalues are sorted lexicographically by real part and then imaginary part.
The first :math:`D/2` columns are anchored at the layer top and the remaining
:math:`D/2` at the layer bottom:

.. math::

   a_{\ell,e}&=\begin{cases}
      t_\ell,&e\le D/2,\\
      t_{\ell+1},&e>D/2,
   \end{cases}\\
   H^{am\ell}_{q e}(\tau)
      &=G^{am\ell}_{q e}\exp[K^{am\ell}_e(\tau-a_{\ell,e})].

For the usual paired spectrum this anchors negative-real-part modes at the top
and positive-real-part modes at the bottom.  The code classifies by sorted
position rather than explicitly checking the sign.  When comparing with a
formulation that writes attenuation as :math:`\exp(-k\tau)`, remember that the
C++ propagation constant has the opposite sign.

The full eigenproblem is intentionally retained instead of using the optional
reduced polarized problem.  It is simpler to audit against the published
equations and supports general Mueller coupling, but it makes VDISORT
substantially more expensive than scalar DISORT.

VDISORT particular-solution right-hand sides
=============================================

For :math:`\epsilon_0=1` and :math:`\epsilon_m=2` when :math:`m>0`, the
direct-beam right-hand side is

.. math::

   b^{am\ell}_{\rho(i,s_o)}
      =\frac{\epsilon_m\omega_\ell}{4\pi\mu_i}
       \left[\boldsymbol P^{am\ell}_{i,\mathrm b}
             \boldsymbol S_{\mathrm b}\right]_{s_o}.

The direct beam travels in direction :math:`-\mu_0`.  Substitution of
:math:`\boldsymbol B^{am\ell}\exp(-\tau/\mu_0)` into the ODE gives the matrix
solve

.. math::

   (\boldsymbol A^{am\ell}+\mu_0^{-1}\boldsymbol I_D)
      \boldsymbol B^{am\ell}=\boldsymbol b^{am\ell}.

For the internal source, let its layer coordinate and polynomial be

.. math::

   x_\ell(\tau)=a_\ell\tau+b_\ell,\qquad
   \boldsymbol q_\ell(x)=\sum_{p=0}^{P-1}\boldsymbol q_{\ell,p}x^p.

Only :math:`m=0` receives this source.  The cosine combined system receives its
:math:`I,Q` components and the sine combined system receives its :math:`U,V`
components.  With

.. math::

   \boldsymbol v_\ell(x)=\sum_{p=0}^{P-1}\boldsymbol V_{\ell,p}x^p,

the coefficient vectors are obtained backwards from
:math:`\boldsymbol V_{\ell,P}=0` using

.. math::

   \boldsymbol A^{a0\ell}\boldsymbol V_{\ell,p}
      =\widehat{\boldsymbol M}^{-1}\boldsymbol q^a_{\ell,p}
       +a_\ell(p+1)\boldsymbol V_{\ell,p+1}.

The factor :math:`a_\ell` comes from differentiating the affine source
coordinate.  The complete particular field is

.. math::

   \boldsymbol p^{am\ell}(\tau)
      =\boldsymbol B^{am\ell}\exp(-\tau/\mu_0)
       +\delta_{m0}\sum_p\boldsymbol V_{\ell,p}x_\ell(\tau)^p.

Scalar limit of the VDISORT equations
=====================================

The VDISORT transfer equations contain the scalar DISORT equation as an
invariant subproblem when all of the following hold:

* radiation is embedded as :math:`[I,0,0,0]^{\mathsf T}`;
* only the :math:`(0,0)` Mueller element contains the scalar phase function;
* sources and boundaries contain only an :math:`I` component; and
* surface reflection maps only incident :math:`I` to outgoing :math:`I`.

Then :math:`Q`, :math:`U`, and :math:`V` remain zero and the :math:`I`
component obeys the scalar equation.  This is only a mathematical and physical
reduction.  The VDISORT implementation does not detect this case, dispatch to
``disort.cpp``, or reduce the dimensions of its matrices.  It still assembles
and diagonalizes the full :math:`4N_q\times4N_q` vector system for both combined
Fourier systems and solves the corresponding polarized boundary problem.
Consequently the scalar embedding is useful for validation, but retains
VDISORT's greater memory and computational cost.

This scalar limit is tested against the same supported numerical reference set
as the scalar solver.

Sources and boundary conditions
*******************************

Both solvers support prescribed upper and lower diffuse boundary fields, a
direct beam, an internal source polynomial in every layer, and a reflecting
lower boundary.  Only the zeroth Fourier mode contains an azimuth-independent
internal source.

In both cores the polynomial represents a source radiance.  VDISORT takes the
Stokes source function :math:`\boldsymbol B(\tau)`, and the transport code
applies :math:`1-\omega`.  Thus an unpolarized thermal source is supplied as

.. math::

   \boldsymbol B(\tau)=[B(\tau),0,0,0]^{\mathsf T},

and VDISORT forms the emission vector
:math:`\boldsymbol q=(1-\omega)\boldsymbol B` internally.  VDISORT can evaluate
its polynomial in a layer-local affine coordinate, which avoids refitting
coefficients when the natural source coordinate is not the global optical
depth.  In both implementations, ``source_poly_coeffs`` retains the physical
source-function input while ``scaled_source_poly_coeffs`` caches the
transport-equation emission coefficients used by the solver.

The direct beam is stored separately from the diffuse quadrature field.  In
VDISORT it is a full Stokes vector and is considered present when its intensity
component is positive.  Its scattering phase matrices are evaluated at the
fixed incident direction :math:`-\mu_0`.

Surface reflection
******************

The lower boundary uses Fourier coefficients of a bidirectional reflectance
distribution function (BRDF).  For a scalar raw BRDF
:math:`\rho(\mu,\mu',\Delta\phi)` in inverse steradians, the helper used by the
core defines

.. math::

   R_m(\mu,\mu') = (2-\delta_{m0})
      \int_0^\pi \rho(\mu,\mu',\varphi)\cos(m\varphi)\,\mathrm d\varphi.

Consequently a Lambertian BRDF :math:`A/\pi` produces :math:`R_0=A`.  The
azimuth integral is evaluated by Gauss--Legendre quadrature.  Hapke, Cox--Munk,
RPV, and Ross--Li raw scalar models use this same transformation.

Scalar DISORT applies the corresponding diffuse-reflection quadrature and
direct-beam reflection in the lower boundary equation.  VDISORT generalizes
each Fourier coefficient to a :math:`4\times4` Mueller operator and retains a
cosine and sine operator for each mode.  The core applies the angular
quadrature weights exactly once; callbacks must therefore return the Fourier
BRDF or Mueller-BRDF coefficient itself, without a quadrature weight folded
into it.

For positive outgoing and incident quadrature streams, VDISORT assembles the
actual bottom reflection matrix as

.. math::

   \left[\boldsymbol{\mathcal R}^{am}
   \right]_{\rho(i,s_o),\rho(j,s_i)}
   =\pi\gamma_m w_j\mu_j
     \left[\boldsymbol R^{\mathrm{raw},am}_{ij}\right]_{s_os_i},
   \qquad
   \gamma_m=\begin{cases}1,&m=0,\\[2pt]\tfrac12,&m>0.\end{cases}

The separately assembled reflected direct-beam vector is

.. math::

   \boldsymbol d^{am}_{\rho(i,s_o)}
   =\frac{\mu_0}{2}\exp(-t_L/\mu_0)
     \left[\boldsymbol R^{\mathrm{raw},am}_{i,\mathrm b}
           \boldsymbol S_{\mathrm b}\right]_{s_o}.

These are the :math:`\boldsymbol{\mathcal R}_m` and
:math:`\boldsymbol d_m` objects used in the global bottom row below.

The built-in raw physical models are scalar.  Embedding them in VDISORT tests
therefore exercises only the Mueller :math:`(0,0)` element.  A genuinely
polarizing surface requires Mueller-valued cosine and sine Fourier modes.

Global boundary-value matrix
****************************

.. important::

   There are two different matrix problems.  The layer eigenproblem
   :math:`\boldsymbol L_{\mathrm{eig}}\boldsymbol g
   =k\boldsymbol R_{\mathrm{eig}}\boldsymbol g` determines
   :math:`\boldsymbol K` and :math:`\boldsymbol G`.  The ``LHS`` and ``RHS``
   objects in ``solve_for_coefs()`` belong to a subsequent global boundary
   solve.  They determine the modal constants, not the eigenvalues.

The scalar problem below has state dimension :math:`d=N_q`; the vector problem
has :math:`d=4N_q`.  In either case let :math:`h=d/2`, let
:math:`\boldsymbol S_+` select the upward half of a state, and let
:math:`\boldsymbol S_-` select the downward half.  For each Fourier system,
write the field in layer :math:`\ell` as

.. math::

   \boldsymbol u_\ell(\tau)
      =\boldsymbol H_\ell(\tau)\boldsymbol c_\ell
       +\boldsymbol p_\ell(\tau),

where :math:`\boldsymbol H_\ell` is the anchored eigenvector/exponential matrix,
:math:`\boldsymbol c_\ell` contains :math:`d` unknown modal constants, and
:math:`\boldsymbol p_\ell` is the sum of all particular solutions.  If
:math:`\boldsymbol{\mathcal R}_m` is the already quadrature-weighted bottom
reflection operator, the unknown and equation ordering is

.. math::

   \boldsymbol c
      =(\boldsymbol c_0,\boldsymbol c_1,\ldots,
        \boldsymbol c_{L-1})^{\mathsf T},

.. math::

   \underbrace{\begin{pmatrix}
     \boldsymbol S_-\boldsymbol H_0(t_0)
       &0&\cdots&0\\
     \boldsymbol H_0(t_1)
       &-\boldsymbol H_1(t_1)&\cdots&0\\
     0&\boldsymbol H_1(t_2)
       &\ddots&\vdots\\
     \vdots&\ddots&\ddots
       &-\boldsymbol H_{L-1}(t_{L-1})\\
     0&\cdots&0&
       (\boldsymbol S_+-\boldsymbol{\mathcal R}_m\boldsymbol S_-)
       \boldsymbol H_{L-1}(t_L)
   \end{pmatrix}}_{\boldsymbol L_{\mathrm{BC}}}
   \boldsymbol c
   =
   \underbrace{\begin{pmatrix}
     \boldsymbol b^\downarrow_m
       -\boldsymbol S_-\boldsymbol p_0(t_0)\\
     \boldsymbol p_1(t_1)-\boldsymbol p_0(t_1)\\
     \vdots\\
     \boldsymbol p_{L-1}(t_{L-1})
       -\boldsymbol p_{L-2}(t_{L-1})\\
     \boldsymbol b^\uparrow_m+\boldsymbol d_m
       -\boldsymbol S_+\boldsymbol p_{L-1}(t_L)
       +\boldsymbol{\mathcal R}_m\boldsymbol S_-
          \boldsymbol p_{L-1}(t_L)
   \end{pmatrix}}_{\boldsymbol R_{\mathrm{BC}}}.

The first :math:`h` rows impose the top downward boundary.  Each internal
interface contributes :math:`d` continuity rows, and the last :math:`h` rows
impose the reflecting bottom boundary.  The count is therefore

.. math::

   h+(L-1)d+h=Ld,

matching the :math:`Ld` unknown modal constants.  Only adjacent layer blocks
couple, so the matrix is stored as a band matrix.  Both implementations use
upper and lower half-bandwidth :math:`3h-1`.  In VDISORT this entire system is
complex and is assembled independently for every :math:`(a,m)`.

The vector :math:`\boldsymbol d_m` is direct-beam reflection.  With no
reflecting surface, :math:`\boldsymbol{\mathcal R}_m=0` and
:math:`\boldsymbol d_m=0`.  The formula also shows why a particular solution
appears as a difference across an interface even though the total radiance is
continuous.

Implementation array layouts
============================

The similarly named arrays in the two cores do not always store the same
factorization:

.. list-table::
   :header-rows: 1
   :widths: 24 35 41

   * - Mathematical object
     - CPP-DISORT
     - VDISORT
   * - Layer transport matrix :math:`\boldsymbol A`
     - Implicit in ``D_pos``, ``D_neg``, ``apb``, and ``amb``; only the
       :math:`N\times N` reduced product ``sqr`` is diagonalized
     - Local :math:`4N_q\times4N_q` real ``A_real`` is copied to a complex
       matrix and diagonalized
   * - Propagation constants :math:`\boldsymbol K`
     - ``K_collect[m, layer, eigen]`` with negative half followed by its
       positive partners
     - ``K_collect[alpha, m, layer, eigen]`` after complex lexicographic sort
   * - Eigenvectors :math:`\boldsymbol G`
     - ``G_collect[m, layer, state, eigen]``
     - ``G_collect[alpha, m, layer, state, eigen]``
   * - Modal constants :math:`\boldsymbol c_\ell`
     - The band-solve ``RHS`` is overwritten by the constants
     - The complex band-solve ``rhs`` is copied into ``GC_collect``
   * - ``GC_collect``
     - Stores :math:`G_{qe}c_e` for every state and eigenmode
     - Despite the name, stores only :math:`c_e`; multiplication by
       :math:`G_{qe}` occurs during field reconstruction
   * - Beam particular solution :math:`\boldsymbol B`
     - ``B_collect[m, layer, stream]``
     - ``B_collect[alpha, m, layer, stream]`` of Stokes vectors
   * - Polynomial particular solution
     - Evaluated through scalar scratch arrays ``SRC0``, ``SRC1``, and
       ``SRCB``
     - Coefficients are Stokes vectors in
       ``source_collect[alpha, m, layer, stream, power]``
   * - Layer-bottom total field
     - ``um[layer, m, stream]``
     - ``um[layer, alpha, m, stream]`` of Stokes vectors

In particular, ``um`` is a cached total quadrature field, not an eigenvector or
a modal-coefficient array.

Delta-M scaling
***************

For a forward-peak fraction :math:`f`, original single-scattering albedo
:math:`\omega`, and normalized removed-peak moments :math:`r_l`, scalar DISORT
uses

.. math::

   s &= 1-\omega f,\\
   \Delta\tau' &= s\,\Delta\tau,\\
   \omega' &= \frac{\omega(1-f)}{1-\omega f},\\
   \chi'_l &= \frac{\chi_l-f r_l}{1-f}.

Classical delta-M has :math:`r_l=1`.  Physical optical depths remain available
to callers, while the eigensystem and direct-beam particular solution use the
scaled coordinate :math:`\tau'`.

Internal-source polynomials retain their physical meaning under this coordinate
change.  In layer :math:`l`, write the cumulative transformation as

.. math::

   \tau' = s_l\tau+c_l,
   \qquad c_l=t_l'-s_l t_l,

where :math:`t_l` and :math:`t_l'` are the physical and scaled optical depths
at the layer top.  The solver composes each physical polynomial with the
inverse affine map once,

.. math::

   B_l'(\tau')=B_l\!\left(\frac{\tau'-c_l}{s_l}\right),

and includes the absorption factor in the derived coefficients.  Since

.. math::

   1-\omega'=\frac{1-\omega}{s_l},

the resulting ``scaled_source_poly_coeffs`` represent the required transformed
emission :math:`(1-\omega)B_l(\tau)/s_l` in the :math:`\tau'` equation.  The
cached coefficients are reused directly by quadrature, flux, gridded, and
arbitrary-angle evaluations.

The delta-M-plus option represents the removed peak by Gaussian Legendre
moments

.. math::

   r_l=\exp\!\left(-\frac{l^2}{2\sigma^2}\right).

The width and fraction are inferred from the first two phase moments beyond
the retained transport expansion.  If any layer fails the consistency checks
for a useful Gaussian tail, all layers fall back to classical delta-M.

VDISORT's transport object receives already transformed Mueller phase
operators and optical depths; it does not invent a general Mueller-valued
delta-M-plus transformation.  The validated delta-M-plus path for VDISORT is
the scalar :math:`M_{00}` reduction.  General polarized delta-M corrections
require an explicitly defined removed Mueller peak in a common laboratory
reference frame.

IMS and TMS intensity corrections
=================================

Delta-M improves the diffuse solution but removes structure from the direct
beam aureole.  The Nakajima--Tanaka correction used here is the sum of a
truncated-multiple-scattering (TMS) term and an improved-multiple-scattering
(IMS) term.  TMS formally integrates the difference between the original and
transport phase functions along the direct-beam path.  Stable limiting kernels
are used where :math:`\mu` approaches :math:`\mu_0`, avoiding cancellation in
expressions containing :math:`1/\mu-1/\mu_0`.

By default, IMS is applied only to downward directions within 10 degrees of
the incident beam and is subtracted from the TMS-corrected field.  The scalar
core also exposes the alternative sign/domain convention for explicit callers.

Scalar IMS/TMS is available at quadrature and arbitrary user angles and has a
gridded path which reuses angle-dependent work.  It is restricted to classical
delta-M.  It is intentionally disabled for a non-classical removed peak,
including delta-M-plus.

VDISORT provides a fully Mueller-valued correction cache.  Its TMS operator is
formed from the difference between original and transport Mueller phase
matrices.  For normalized removed peak :math:`\mathbf R`, the IMS angular
operator contains

.. math::

   2\mathbf R(\Omega,\Omega_0)
   - \frac{1}{4\pi}\int_{4\pi}
       \mathbf R(\Omega,\Omega')\mathbf R(\Omega',\Omega_0)
       \,\mathrm d\Omega'.

The matrix product preserves polarization and reference-frame rotations.  This
angular convolution is expensive, so it is computed once for a selected set of
user directions and azimuths and cached.  A supplied analytic pair-convolution
can replace the numerical intermediate-angle quadrature.  Evaluation at
different optical depths then reuses the cached operators.

Arbitrary-angle radiances
*************************

Quadrature-stream radiances come directly from the discrete-ordinate field.
Radiances at other nonzero direction cosines are obtained from the formal
solution along the requested ray; they are not merely interpolated values of
the final radiance.

Scalar DISORT reconstructs the angular source from the quadrature solution.
It uses barycentric interpolation for the required angular terms and analytic,
cancellation-safe exponential integrals through every crossed layer.

For polarized radiation, a phase matrix sampled only at the quadrature output
directions is insufficient to reconstruct a general outgoing direction.
VDISORT therefore requires the directional diffuse and direct-beam Mueller
phase matrices at every requested user direction.  The supplied matrices are
contracted with the stored eigenvectors, integration constants, beam
particular solution, and polynomial particular solution.  The resulting
layerwise complex exponentials and exponential--polynomial terms are integrated
analytically with cancellation-safe limiting functions, following the same
formal-solution structure as scalar DISORT.  A bulk call forms the
azimuth-independent combined Fourier modes locally for all requested optical
depths, then synthesizes every requested azimuth; no persistent user-radiance
cache is part of the solver state.

The diffuse top and bottom boundary fields are stored only at quadrature
directions.  Their values at a new direction are therefore not uniquely
defined by the core inputs.  The test-facing arbitrary-angle path uses stable
half-range barycentric interpolation, which is exact for the constant boundary
fields exercised by the reference suite.  General polarized user-angle use
would require explicit directional boundary data or a boundary callback.

The direction :math:`\mu=0` is excluded in both formal solutions because the
ray equation contains :math:`1/\mu`.

Fluxes and heating-rate integrand
*********************************

Only the zeroth Fourier mode contributes to hemispheric flux.  Both cores
return positive upward, positive diffuse-downward, and positive direct-downward
fluxes.  With the positive quadrature nodes,

.. math::

   F^\uparrow &= 2\pi\sum_{i=1}^{N}w_i\mu_i I^+_i,\\
   F^\downarrow_{\mathrm{diff}} &=
      2\pi\sum_{i=1}^{N}w_i\mu_i I^-_i,\\
   F^\downarrow_{\mathrm{dir}} &=
      \mu_0 I_{\mathrm b}\exp(-\tau/\mu_0).

In VDISORT these quantities use the :math:`I` component of each Stokes vector;
the core does not currently return hemispherically integrated :math:`Q`,
:math:`U`, or :math:`V` fluxes.

The combined flux call also returns the local flux-divergence or heating-rate
integrand called ``DFDT`` by DISORT.  If

.. math::

   J = \frac{1}{4\pi}\int_{4\pi}I\,\mathrm d\Omega

includes the direct beam, scalar DISORT evaluates

.. math::

   \mathrm{DFDT}=4\pi(1-\omega)(J-B).

VDISORT uses the same expression because its polynomial stores
:math:`\boldsymbol B`, not the already absorption-weighted emission vector.
Computing all flux quantities in one call reuses a single zeroth-mode radiance
evaluation.

Numerical behavior and limitations
**********************************

The following details are intentional and should be considered when changing
or comparing the cores:

* **Conservative scattering.**  Exact :math:`\omega=1` gives a zero eigenvalue
  in the scalar zeroth-mode reduction.  Scalar DISORT internally replaces it
  by :math:`1-10^{-8}`, the smallest tested dither which retains the original
  DISORT reference results with this eigensolver.  VDISORT accepts
  :math:`\omega=1` directly, but conservative vector eigensystems can be poorly
  conditioned.
* **Complex VDISORT modes.**  A real physical solution can use complex
  eigenpairs.  The reconstructed imaginary part is required to cancel to a
  relative tolerance of about :math:`2\,10^{-8}` before the real part is
  returned.  Failure indicates an ill-conditioned or incorrectly assembled
  system, not a physical complex radiance.
* **VDISORT spectral split.**  Anchoring assumes that sorting by real part
  places exactly half the modes on each side of the spectrum.  Degenerate or
  nearly imaginary modes can make this classification ill-conditioned because
  the implementation splits by sorted position, not by an explicit sign test.
* **Beam--eigenmode resonance.**  The direct-beam particular solution contains
  :math:`(K_e+1/\mu_0)^{-1}`.  If a propagation constant is close to
  :math:`-1/\mu_0`, the direct-beam solve is poorly conditioned.  The cores do
  not currently issue a proximity warning.
* **Transparent layers.**  A layer with exactly zero optical thickness and no
  interaction can produce a degenerate eigensystem.  Reference surface-only
  tests use a negligible but nonzero optical thickness and albedo instead of
  claiming exact transparent-layer support.
* **No absorption cutoff.**  An absorption-optical-depth shortcut is not
  implemented.  The complete physical solution is evaluated instead.
* **No special-boundary shortcut.**  ``IBCND=1`` albedo/transmission is not a
  separate algorithm here.  The same quantities are obtained from the regular
  boundary-value solution.
* **No separate two-stream solver.**  Setting :math:`N_q=2` selects a
  two-stream discrete-ordinate instance rather than a distinct approximate
  two-stream algorithm.
* **Fourier truncation.**  Features narrow in relative azimuth require enough
  atmospheric and surface Fourier modes.  Likewise, sharply peaked phase
  functions require enough Legendre moments or a suitable delta-M treatment.
* **Phase conventions.**  VDISORT assumes all Mueller matrices have consistent
  propagation directions and Stokes reference planes.  The supplied combined
  matrices already include the cosine/sine transformation.  Applying that
  transformation twice, adding an extra factor of two, or mixing
  :math:`[I,Q,U,V]` with :math:`[I_v,I_h,U,V]` changes the physics.
* **Update cost.**  Constructing or updating either solver performs the
  eigendecompositions and global boundary solves.  Evaluation at cached
  quadrature points is much cheaper.  VDISORT intentionally contains no
  OpenMP parallel region in these core algorithms; ARTS normally parallelizes
  independent frequencies outside the solver.
* **Thread use.**  A fully constructed solver may be shared for read-only
  evaluation when each caller owns its scratch objects.  Updating the solver
  or sharing scratch storage concurrently is unsafe.

Validation scope
****************

The scalar core is checked against published numerical reference values for
the supported plane-parallel test problems, including arbitrary angles, fluxes,
DFDT, delta-M corrections, physical BRDFs, conservative scattering, multilayer
continuity, and delta-M-plus.  The pseudo-spherical test case is intentionally
unsupported.

The strict scalar embedding of VDISORT is checked against the same supported
reference problems, with :math:`Q=U=V=0` asserted.  These tests establish the
normalization and scalar limit of the vector equations.  Analytic polarized
two-stream tests additionally cover :math:`I/Q` coupling, complex :math:`U/V`
eigenpairs, a polarized direct beam, polarized absorption, vector internal
sources, and a polarized reflecting boundary.  These focused tests do not
replace comparison with an independent general-purpose polarized reference,
particularly for reference-plane rotations, many-stream Mueller problems, and
genuinely polarized IMS/TMS corrections.
