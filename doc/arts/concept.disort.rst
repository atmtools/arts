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
<https://doi.org/10.1364/AO.27.002502>`_.  The scalar implementation is
structurally based on `PythonicDISORT
<https://doi.org/10.21105/joss.06442>`_, but has been extended to reproduce the
supported DISORT 4.0.99 reference problems.
The vector implementation follows the full-eigenproblem formulation of
`Lin et al. (2022) <https://doi.org/10.3389/frsen.2022.880768>`_ while retaining
the organization of the scalar solver.

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

The geometry is deliberately plane-parallel.  Neither core implements the
pseudo-spherical direct-beam approximation from DISORT Problem 16.  Curved ray
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
   = \sum_{l=0}^{L}(2l+1)\chi_l P_l(\cos\Theta),
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

The scalar code exploits the symmetry between positive and negative streams.
It reduces the :math:`N_q`-dimensional eigensystem to an
:math:`N_q/2`-dimensional real eigensystem, takes the paired propagation
constants :math:`\pm k`, and reconstructs the full eigenvectors.  Homogeneous,
direct-beam, and polynomial particular solutions are then joined by continuity
at every layer interface.  The upper and lower boundary conditions form a
banded linear system covering all layers.

The exponentials belonging to modes which grow with increasing optical depth
are anchored at the lower edge of a layer; the others are anchored at its
upper edge.  This avoids explicitly evaluating the largest growing
exponentials across a complete layer.

Polarized radiative transfer
****************************

VDISORT transports the ARTS Stokes vector

.. math::

   \boldsymbol I = [I,Q,U,V]^{\mathsf T},

not the alternative :math:`[I_v,I_h,U,V]^{\mathsf T}` representation.  The
vector equation is

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
cosine and sine Mueller coefficients according to Eqs. 78--82 of Lin et al.
The physical field is reconstructed schematically as

.. math::

   \begin{aligned}
   [I,Q]^{\mathsf T}
     &= \sum_m\left(\boldsymbol C^m_{IQ}\cos m\Delta\phi
                    +\boldsymbol S^m_{IQ}\sin m\Delta\phi\right),\\
   [U,V]^{\mathsf T}
     &= \sum_m\left(\boldsymbol S^m_{UV}\cos m\Delta\phi
                    +\boldsymbol C^m_{UV}\sin m\Delta\phi\right),
   \end{aligned}

with :math:`\Delta\phi=\phi_0-\phi`.  The swapped roles in the second line are
important; treating every Stokes component like scalar intensity gives the
wrong signs and azimuthal behavior.

Scattering data produced in a complex spectral representation are split as

.. math::

   \mathbf P_{\cos}=\Re\mathbf P_{\mathbb C},\qquad
   \mathbf P_{\sin}=-\Im\mathbf P_{\mathbb C}.

No factor of two is introduced by this split.  Fourier normalization is handled
when the modes are combined and synthesized.  This interpretation applies to
the phase-matrix output, rather than changing the scattering-species data
model.

For each layer, Fourier order, and combined system, VDISORT solves the full
:math:`4N_q\times4N_q` real, generally non-symmetric transport eigenproblem.
Its eigenvalues and eigenvectors may be complex even though the transport
matrix and physical radiance are real.  Conjugate pairs must cancel during
reconstruction.  Eigenmodes are sorted and half are anchored at each layer
edge before the global complex banded boundary system is solved.

The full eigenproblem is intentionally retained instead of using the optional
reduced polarized eigenproblem.  It is simpler to audit against the published
equations and supports general Mueller coupling, but it makes VDISORT
substantially more expensive than scalar DISORT.

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

This scalar limit is tested against the same supported DISORT 4.0.99 reference
problems as the scalar solver.

Sources and boundary conditions
*******************************

Both solvers support prescribed upper and lower diffuse boundary fields, a
direct beam, an internal source polynomial in every layer, and a reflecting
lower boundary.  Only the zeroth Fourier mode contains an azimuth-independent
internal source.

There is one notable normalization difference inside the cores.  Scalar
DISORT's polynomial represents a source radiance :math:`B(\tau)` and the
transport code applies :math:`1-\omega`.  VDISORT's polynomial is the complete
emission vector :math:`\boldsymbol q(\tau)` already appearing in the transfer
equation.  Thus an unpolarized thermal source is embedded as

.. math::

   \boldsymbol q(\tau)
   = (1-\omega)[B(\tau),0,0,0]^{\mathsf T}.

This distinction must be observed by adapters comparing the two solvers.
VDISORT can evaluate its polynomial in a layer-local affine coordinate, which
avoids refitting coefficients when the natural source coordinate is not the
global optical depth.

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

The built-in raw physical models are scalar.  Embedding them in VDISORT tests
therefore exercises only the Mueller :math:`(0,0)` element.  A genuinely
polarizing surface requires Mueller-valued cosine and sine Fourier modes.

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

The delta-M-plus option represents the removed peak by Gaussian Legendre
moments

.. math::

   r_l=\exp\!\left(-\frac{l^2}{2\sigma^2}\right).

The width and fraction are inferred from the first two phase moments beyond
the retained transport expansion, following DISORT 4.0.99.  If any layer fails
the DISORT guards for a useful Gaussian tail, all layers fall back to classical
delta-M, also matching the reference implementation.

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

Under the DISORT convention, IMS is applied only to downward directions within
10 degrees of the incident beam and is subtracted from the TMS-corrected field.
The scalar core also retains the alternative PythonicDISORT sign/domain
convention for explicit callers, but DISORT reference comparisons use the
DISORT convention.

Scalar IMS/TMS is available at quadrature and arbitrary user angles and has a
gridded path which reuses angle-dependent work.  It is restricted to classical
delta-M.  It is intentionally disabled for a non-classical removed peak,
including delta-M-plus, matching DISORT 4.0.99.

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
cancellation-safe exponential integrals through every crossed layer.  This is
the counterpart of the ``TERPEV``, ``TERPSO``, and ``USRINT`` path in original
DISORT.

For polarized radiation, a phase matrix sampled only at the quadrature output
directions is insufficient to reconstruct a general outgoing direction.
VDISORT therefore requires the directional diffuse and direct-beam Mueller
phase matrices at every requested user direction.  It formally integrates the
source with segmented 12-point Gauss--Legendre quadrature.  The resulting
azimuth-independent combined Fourier modes can be cached and synthesized at
many azimuths without repeating the depth integration.

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

VDISORT uses the equivalent form

.. math::

   \mathrm{DFDT}=4\pi\left[(1-\omega)J-q_I\right],

because its polynomial already stores the complete emission term
:math:`q_I=(1-\omega)B`.  Computing all flux quantities in one call reuses a
single zeroth-mode radiance evaluation.

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
* **Transparent layers.**  A layer with exactly zero optical thickness and no
  interaction can produce a degenerate eigensystem.  Reference surface-only
  tests use a negligible but nonzero optical thickness and albedo instead of
  claiming exact transparent-layer support.
* **No absorption cutoff.**  The original DISORT absorption-optical-depth
  shortcut is not implemented.  The complete physical solution is evaluated
  instead.
* **No special-boundary shortcut.**  ``IBCND=1`` albedo/transmission is not a
  separate algorithm here.  The same quantities are obtained from the regular
  boundary-value solution.
* **No separate two-stream solver.**  Setting :math:`N_q=2` selects a
  two-stream discrete-ordinate instance.  It is not the distinct ``TWOSTR``
  approximation found in cDISORT.
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

The scalar core is checked against the original DISORT 4.0.99 numerical
references for Problems 1--15 and 17, including arbitrary angles, fluxes,
DFDT, delta-M corrections, physical BRDFs, conservative scattering, multilayer
continuity, and delta-M-plus.  Problem 16 is the intentionally unsupported
pseudo-spherical case.

The strict scalar embedding of VDISORT is checked against the same supported
reference problems, with :math:`Q=U=V=0` asserted.  These tests establish the
normalization and scalar limit of the vector equations.  Analytic polarized
two-stream tests additionally cover :math:`I/Q` coupling, complex :math:`U/V`
eigenpairs, a polarized direct beam, polarized absorption, vector internal
sources, and a polarized reflecting boundary.  These focused tests do not
replace comparison with an independent general-purpose polarized reference,
particularly for reference-plane rotations, many-stream Mueller problems, and
genuinely polarized IMS/TMS corrections.
