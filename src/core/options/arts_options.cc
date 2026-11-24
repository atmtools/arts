#include "arts_options.h"

#include <auto_species_enum_info.h>
#include <enums-common-helper.h>
#include <mystring.h>

#include <algorithm>
#include <array>
#include <format>
#include <limits>
#include <ranges>
#include <set>
#include <utility>
#include <vector>

using namespace std::literals;

using Value = std::vector<std::string>;

namespace stdr = std::ranges;
namespace stdv = std::views;

namespace {
void fix_static(std::vector<EnumeratedOption>& opts) {
  for (auto& opt : opts) {
    if (opt.desc.empty()) throw std::runtime_error("No empty docstrings");

    if (opt.desc.back() != '\n') opt.desc.push_back('\n');

    for (auto& v : opt.values_and_desc) {
      if (v.size() < 2) throw std::runtime_error("No empty value descriptions");

      auto& l = v.back();

      if (l.empty()) throw std::runtime_error("No empty value descriptions");

      if (l.back() != '\n') v.back().push_back('\n');
    }
  }
}

std::vector<EnumeratedOption> internal_options_create() {
  std::vector<EnumeratedOption> opts;

  opts.emplace_back(EnumeratedOption{
      .name = "OptimalEstimationMethod",
      .desc =
          R"--(Algorithm and linear solver used by optimal estimation.

Every option minimizes the same cost function, so the choice is about how the
minimum is approached and not about what is being minimized:

.. math::
  \chi^2 = \chi^2_y + \chi^2_x,

.. math::
  \chi^2_y = \frac{1}{m}\left(\vec{y}-\vec{y}_f\right)^\top
             \mathbf{S}_\epsilon^{-1}\left(\vec{y}-\vec{y}_f\right),
  \qquad
  \chi^2_x = \frac{1}{m}\left(\vec{x}-\vec{x}_a\right)^\top
             \mathbf{S}_a^{-1}\left(\vec{x}-\vec{x}_a\right).

The a priori term belongs to the cost for every option, and no option changes
the statistical meaning of the two covariance matrices.  What an option selects
is how the state is stepped and which linear system produces that step.  Both
are built from the half-gradient and the Gauss--Newton approximation of the
half-Hessian of :math:`m\chi^2`:

.. math::
  \vec{g} = \mathbf{J}^\top\mathbf{S}_\epsilon^{-1}
            \left(F\left(\vec{x}\right)-\vec{y}\right)
            + \mathbf{S}_a^{-1}\left(\vec{x}-\vec{x}_a\right),
  \qquad
  \mathbf{H} = \mathbf{J}^\top\mathbf{S}_\epsilon^{-1}\mathbf{J}
               + \mathbf{S}_a^{-1}.

There are two ways to step.  ``li`` and ``gn`` take the undamped
Gauss--Newton step

.. math::
  \mathbf{H}\Delta\vec{x} = -\vec{g},

which ``li`` takes exactly once and ``gn`` repeats.  ``lm`` instead takes the
damped step of :ref:`sec-oem-damping`

.. math::
  \left(\mathbf{H}+\gamma\mathbf{D}\right)\Delta\vec{x} = -\vec{g},
  \qquad
  \mathbf{D} = \operatorname{diag}
               \left(\operatorname{diag}\left(\mathbf{S}_a^{-1}\right)\right),

where :math:`\gamma` is raised on a rejected trial and lowered on an accepted
one.  At :math:`\gamma=0` this is the Gauss--Newton step again.

There are also two spaces to solve in.  The plain options solve the
:math:`n \times n` state-space system written above.  The ``_m`` options solve
the :math:`m \times m` measurement-space system instead:

.. math::
  \mathbf{M} = \mathbf{J}\mathbf{S}_a\mathbf{J}^\top + \mathbf{S}_\epsilon,
  \qquad
  \vec{r} = \vec{y} - F\left(\vec{x}_i\right)
            + \mathbf{J}\left(\vec{x}_i-\vec{x}_a\right),

.. math::
  \mathbf{M}\vec{u} = \vec{r},
  \qquad
  \vec{x}_{i+1} = \vec{x}_a + \mathbf{S}_a\mathbf{J}^\top\vec{u}.

The two systems have the same solution in exact arithmetic, so the choice is
one of size and conditioning: the state-space system is the smaller one when
there are more measurements than state elements, and the measurement-space
system when there are fewer.  Finally, the ``_cg`` options replace the direct
solve of the selected system with a conjugate-gradient solve, which never
forms the inverse and is worth its tolerance only when the direct solve is the
bottleneck.

where:

- :math:`\vec{y}` is the measurement vector,
- :math:`\vec{y}_f = F\left(\vec{x}\right)` is the simulated measurement,
- :math:`\vec{x}` is the model state and :math:`\vec{x}_a` its a priori,
- :math:`F` is the forward model, run through *inversion_iterate_agenda*,
- :math:`\mathbf{J}` is the Jacobian of :math:`F` at the current state,
- :math:`\mathbf{S}_\epsilon` is the measurement error covariance,
- :math:`\mathbf{S}_a` is the a priori covariance of the model state,
- :math:`m` is the number of measurements, and
- :math:`n` is the number of retrieved state elements.

Four limits end a retrieval, and which of them can be reached depends on the
option.  An iteration budget stops the outer loop of the iterated options and
reports ``IterationLimit``; the linear options report it after their single
step.  A convergence limit compares the state-step measures of
:ref:`sec-oem-convergence` against ``stop_dx``, which is the only way to
report ``Converged``.  A conjugate-gradient budget applies to the ``_cg``
options only: an unconverged step is rejected, and ``li`` and ``gn`` then stop
with ``LinearSolverLimit`` while ``lm`` first retries with more damping.  A
damping limit applies to ``lm`` alone, which reports ``DampingLimit`` when it
runs out of damping without an acceptable step.  A starting cost above
``max_start_cost`` skips the inversion for every option and reports
``StartCostLimit``.

.. note::
  ``lm`` has no measurement-space form, as the damping is applied to
  :math:`\mathbf{H}` in state space.

  ``ml`` and ``ml_cg`` are historical spellings of ``lm`` and ``lm_cg``.  They
  select neither maximum-likelihood estimation nor a retrieval without the a
  priori term.
)--",
      .values_and_desc =
          {
              Value{"li", "Linear", "A single undamped step, solving the state-space system directly."},
              Value{"li_m",
                    "LinearMeasurementSpace",
                    "A single undamped step, solving the measurement-space system directly."},
              Value{"li_cg",
                    "LinearConjugateGradient",
                    "A single undamped step, solving the state-space system by conjugate gradient."},
              Value{"li_cg_m",
                    "LinearConjugateGradientMeasurementSpace",
                    "A single undamped step, solving the measurement-space system by conjugate gradient."},
              Value{"gn", "GaussNewton", "Iterated undamped steps, solving the state-space system directly."},
              Value{"gn_m",
                    "GaussNewtonMeasurementSpace",
                    "Iterated undamped steps, solving the measurement-space system directly."},
              Value{"gn_cg",
                    "GaussNewtonConjugateGradient",
                    "Iterated undamped steps, solving the state-space system by conjugate gradient."},
              Value{"gn_cg_m",
                    "GaussNewtonConjugateGradientMeasurementSpace",
                    "Iterated undamped steps, solving the measurement-space system by conjugate gradient."},
              Value{"lm", "ml", "Iterated damped steps, solving the state-space system directly."},
              Value{"lm_cg", "ml_cg", "Iterated damped steps, solving the state-space system by conjugate gradient."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "OptimalEstimationStatus",
      .desc = "Outcome of an optimal-estimation retrieval, reported by OptimalEstimationDiagnostics.\n",
      .values_and_desc =
          {
              Value{"NotRun", "No retrieval has run."},
              Value{"Converged", "The convergence criterion was met or LM established numerical stationarity."},
              Value{"IterationLimit",
                    "The outer iteration budget was reached. Linear methods may report this after their single step."},
              Value{"LinearSolverLimit", "CG exhausted its iteration budget; the unconverged step was rejected."},
              Value{"DampingLimit", "LM reached its damping limit without finding an acceptable step."},
              Value{"Error", "An error was caught during inversion; inspect the diagnostic messages."},
              Value{"StartCostLimit", "The starting cost exceeded max_start_cost and inversion was skipped."},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "AntennaType",
      .desc =
          R"(A switch controlling how monte carlo antenna patterns are handled.
)",
      .values_and_desc = {Value{"PencilBeam", "Pencil beam antenna"},
                          Value{"Gaussian", "Gaussian beam antenna"},
                          Value{"Lookup", "Lookup table antenna"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "RadarUnit",
      .desc            = R"(The unit of a simulated radar return.

The native return is the gate-averaged, attenuated backscatter coefficient.
The equivalent reflectivity factor is

.. math::

    Z_e = \frac{4 \cdot 10^{18} \lambda^4}{\pi^4 |K|^2} \beta,

where :math:`\lambda` is the wavelength [m], :math:`|K|^2` is the reference
dielectric factor, and :math:`\beta` is the native return.
)",
      .values_and_desc = {Value{"unit", "1", "Backscatter coefficient [1 / (m sr)]"},
                          Value{"Ze", "ze", "Equivalent reflectivity factor [mm :math:`^{6}` / m :math:`^{3}`]"},
                          Value{
                              "dBZe", "dbze", "Equivalent reflectivity factor in decibels, :math:`10 \\log_{10} Z_e`"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "RadarRangeMode",
      .desc            = R"(The coordinate of radar range-gate limits.
)",
      .values_and_desc = {Value{"Altitude", "Altitude [m]"},
                          Value{"Distance", "One-way geometric distance from the radar [m]"},
                          Value{"RoundTripTime",
                                "Round-trip time from the radar, including path group refractive indices [s]"},
                          Value{"Legacy",
                                "The ARTS2 rule: Altitude if the largest limit exceeds 1, otherwise RoundTripTime"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "RadarAuxiliary",
      .desc = R"(Auxiliary quantities of a radar simulation.
)",
      .values_and_desc =
          {
              Value{"RadiativeBackground",
                    "Radiative background",
                    "Type of the path end: 0 for space, 1 for the surface, and 2 otherwise"},
              Value{"Backscattering",
                    "backscattering",
                    "Gate-averaged backscatter without attenuation, in the unit of the measurement"},
              Value{"AbsSpeciesExtinction",
                    "Abs species extinction",
                    "Gate-averaged extinction by absorbing species [1 / m]"},
              Value{"ParticleExtinction",
                    "Particle extinction",
                    "Gate-averaged particle extinction, after any extinction scaling [1 / m]"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name = "ray_path_observer_agendaSetGeometricMaxStep",
      .desc =
          R"(For use with *ray_path_observer_agendaSetGeometric*.  Determines how to densify the geometric path.
)",
      .values_and_desc = {Value{"half", "1/2", "Use *ray_pathFillGeometricHalfStep*"},
                          Value{"step", "linear", "Use *ray_pathFillGeometricStepwise*"},
                          Value{"None", "0", "Do not use any step filling method"}},
      .preferred_print = 1,
  });

  opts.emplace_back(EnumeratedOption{
      .name = "settings_agenda_setup_layer_emission_type",
      .desc =
          R"(For atmospheric emission settings with *disort_settings_agendaSetup* and *vdisort_settings_agendaSetup*.
)",
      .values_and_desc = {Value{"None", "No atmospheric emission"},
                          Value{"LinearInTau", "Layer emission is linear in optical depth"},
                          Value{"LinearInTauNonLTE",
                                "Layer emission is linear in optical depth taking non-LTE into account"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "settings_agenda_setup_sun_type",
      .desc =
          R"(For direct radiation settings with *disort_settings_agendaSetup* and *vdisort_settings_agendaSetup*.
)",
      .values_and_desc = {Value{"None", "No direct radiation"},
                          Value{"Sun", "Add a *Sun* object as the direct radiation source"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "settings_agenda_setup_space_type",
      .desc =
          R"(For space boundary settings with *disort_settings_agendaSetup* and *vdisort_settings_agendaSetup*.
)",
      .values_and_desc = {Value{"None", "No isotropic emission from space"},
                          Value{"CosmicMicrowaveBackgroundRadiation",
                                "Isotropic emission as per cosmic microwave background radiation"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "settings_agenda_setup_scattering_type",
      .desc =
          R"(For atmospheric scattering settings with *disort_settings_agendaSetup* and *vdisort_settings_agendaSetup*.
)",
      .values_and_desc = {Value{"None", "No atmospheric scattering"},
                          Value{"ScatteringSpecies", "Use *ArrayOfScatteringSpecies* for scattering"},
                          Value{"ScatteringSpeciesDeltaMPlus",
                                "Use *ArrayOfScatteringSpecies* with spectral delta-M-plus scaling (DISORT only)"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "settings_agenda_setup_surface_type",
      .desc =
          R"(For surface boundary settings with *disort_settings_agendaSetup* and *vdisort_settings_agendaSetup*.
)",
      .values_and_desc = {Value{"None", "No surface emission or scattering"},
                          Value{"Thermal", "Thermal emission, no scattering"},
                          Value{"ThermalLambertian", "Thermal emission, Lambertian scattering"},
                          Value{"Lambertian", "No surface emission, Lambertian scattering"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "SensorKeyType",
      .desc =
          R"(A key for identifying a sensor property
)",
      .values_and_desc = {Value{"freq", "f", "Frequency", "frequency", "Frequency of the sensor"},
                          Value{"zen", "za", "Zenith", "zenith", "Zenith angle of the sensor"},
                          Value{"azi", "aa", "Azimuth", "azimuth", "Azimuth angle of the sensor"},
                          Value{"alt", "h", "Altitude", "altitude", "Altitude of the sensor"},
                          Value{"lat", "la", "Latitude", "latitude", "Latitude of the sensor"},
                          Value{"lon", "lo", "Longitude", "longitude", "Longitude of the sensor"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "HitranLineStrengthOption",
      .desc =
          R"(The way line strength is computed in ARTS when reading Hitran data.

ARTS uses Einstein A-coefficients to compute the line strength. Hitran provides
both the line strength and the Einstein A-coefficient. There is a 1-to-1 conversion
between these two.  However, as with all data, the numbers might differ slightly even
if good cases, so we have provide this selection mechanism to make them match.
)",
      .values_and_desc = {Value{"S", "s", "strength", "s0", "S0", "Line strength"},
                          Value{"A", "a", "einstein", "ein", "A0", "Absorption intensity"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "HitranZeemanOption",
      .desc =
          R"(The way the Zeeman g-factors are set when reading Hitran data.

Reading the g-factors does not turn on the Zeeman effect.
)",
      .values_and_desc = {Value{"None", "No g-factors, they are all 0"},
                          Value{"Compute", "Compute the g-factors, throw if they cannot be computed for a line"},
                          Value{"ComputeOrZero",
                                "Compute the g-factors, set them to 0 if they cannot be computed for a line"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "HitranFileFormatType",
      .desc =
          R"(A file format used to generate an entry in a Hitran-format file.

Note that ARTS do not support all Hitran formats, so if you need a specific
format not listed here, please contribute.

If you download data but do not care for using the parameter in ARTS, use the "unknown" option.

The standard HAPI separator of comma (",") is presumed between fields.
)",
      .values_and_desc =
          {
              Value{"par", "par_line", "The 160-character standard Hitran line format"},
              Value{"statep", "qnsp", "The upper state"},
              Value{"statepp", "qnspp", "The lower state"},
              Value{"unknown", "X", "Skip this field"},
          },
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "SurfaceKey",
      .desc            = R"(A key to identify a surface property.

This is used to identify core components in the surface.
See *SurfaceField* and *SurfacePoint* for its usage.

A core component must either be available or internal calculations will assume
it has a value of 0 once the data is extracted.  This might yield NaN or Inf
values in relevant calculations.  Or it might simply disable a functionality.
)",
      .values_and_desc = {Value{"h", "elevation", "Altitude [m]"}, Value{"t", "temperature", "Temperature [K]"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "AtmKey",
      .desc            = R"(A key to identify an atmospheric property.

This is used to identify core components in the atmosphere.
See *AtmField* and *AtmPoint* for its usage.

A core component must either be available or internal calculations will assume
it has a value of 0 once the data is extracted.  This might yield NaN or Inf
values in relevant calculations.  Or it might simply disable a functionality.
)",
      .values_and_desc = {Value{"t", "temperature", "Temperature [K]"},
                          Value{"p", "pressure", "Pressure [Pa]"},
                          Value{"wind_u", "WindU", "Wind field U-component [m/s]"},
                          Value{"wind_v", "WindV", "Wind field V-component [m/s]"},
                          Value{"wind_w", "WindW", "Wind field W-component [m/s]"},
                          Value{"mag_u", "MagU", "Magnetic field U-component [T]"},
                          Value{"mag_v", "MagV", "Magnetic field V-component [T]"},
                          Value{"mag_w", "MagW", "Magnetic field W-component [T]"},
                          Value{"rot_u", "RotU", "Rotation U-component [m/s], added to wind_u for the Doppler shift"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "SubsurfaceKey",
      .desc            = R"(A key to identify a subsurface property.

This is used to identify core components in the subsurface.
See *SubsurfaceField* and *SubsurfacePoint* for its usage.

A core component must either be available or internal calculations will assume
it has a value of 0 once the data is extracted.  This might yield NaN or Inf
values in relevant calculations.  Or it might simply disable a functionality.
)",
      .values_and_desc = {Value{"t", "temperature", "Temperature [K]"}, Value{"rho", "density", "Density [kg/m^3]"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "MeasurementVectorSumupKernel",
      .desc =
          R"(The kernel to use for summing up the measurement vectors.

This is used when summing up *measurement_vec* and *measurement_jac* from
the observation elements in *measurement_sensors*.  Every kernel calls
*spectral_rad_observer_agenda* once per geometry of each sensor, and gives the
same result up to the rounding of the sums.  The kernels differ in what they
keep in memory and in what they do besides the agenda calls.

The estimates below use, for each sensor in *measurement_sensors*:

- :math:`G`: its geometries, the agenda calls it makes.
- :math:`F`: its frequencies.
- :math:`E`: its elements.
- :math:`W`: its non-zero weights.
- :math:`P`: its pairs of element and geometry where the element has weights.
  At most :math:`W`; for an antenna, the elements times the geometries the
  antenna samples per placement.

The ``report()`` of *measurement_sensors* lists :math:`G`, :math:`F`,
:math:`E`, and :math:`W` of each sensor.  :math:`J` is ``x_size()`` of
*jac_targets*, the number of columns of *measurement_jac*.  :math:`C` is the
number of cores used, at most the total number of geometries.  A number takes
8 bytes, and a Stokes vector 32 bytes.

All kernels keep *measurement_vec* and *measurement_jac*, :math:`8(J+1)` bytes
per element of all sensors.  Each agenda call also uses the working memory of
*spectral_rad_observer_agenda*, on every core that runs one, which comes on top
of the estimates below.  It depends on the agenda, and often dominates.

"Low Memory" keeps, on each core, the spectral radiance and its Jacobian of
the geometry the core simulates: :math:`32(J+1)F` bytes per core.  Besides the
agenda calls, it searches the weights of every element for every geometry,
:math:`GE` searches per sensor, and it sums the weights, :math:`W(J+1)`
products per sensor, one core at a time.  It uses the least memory, but its
searches can dominate when a sensor has both many geometries and many
elements, and the agenda calls are cheap.

"High Performance" keeps the spectral radiance and its Jacobian of every
geometry of every sensor: :math:`32(J+1)GF` bytes per sensor.  It sums the
weights after all agenda calls, in parallel over the elements.  It is the
fastest kernel when this memory is available and the summing is heavy, such as
for a single geometry with wide channels, or for many elements sharing each
geometry.

"Indexed", the default, keeps what "Low Memory" keeps, :math:`32(J+1)F` bytes
per core.  Before the agenda calls, it lists for each geometry the elements
with weights there, and where those weights are.  The list of a sensor takes
about :math:`24P + 8(G+E)` bytes, and building it reads each weight once.  Each
core sums the elements of its geometry after its agenda call.  An element with
weights at that geometry only is written directly.  For each of the other
elements, the core keeps its contribution, :math:`8(J+1)` bytes, and adds it to
*measurement_vec* and *measurement_jac* one core at a time: :math:`P(J+1)`
additions per sensor.

How the kernels grow with the number of cores :math:`C` follows from the
above.  Let :math:`A` be the time of all agenda calls on one core, :math:`B`
the time of the :math:`W(J+1)` products of summing the weights on one core,
and :math:`S` the time of what a kernel does one core at a time.  Sums are over
all sensors, :math:`F` in a per-core term is that of the sensor with the most
frequencies, and :math:`E_g` is the most elements with weights at one
geometry.

- "Low Memory": memory :math:`32(J+1)FC` bytes, linear in :math:`C`.  Its
  :math:`S` is its :math:`GE` searches and :math:`B`, so its time is about
  :math:`\max((A + S)/C, S)`.  It stops gaining from more cores at about
  :math:`C = A/S`.
- "Indexed": memory :math:`(32F + 8E_g)(J+1)C` bytes plus its lists, linear in
  :math:`C` and a little steeper than "Low Memory".  Its :math:`S` is its
  :math:`P(J+1)` additions, so its time is about :math:`\max((A + B + S)/C,
  S)`, plus reading every weight once for the lists.  It stops gaining at about
  :math:`C = (A + B)/S`, later than "Low Memory", as its :math:`S` is smaller.
- "High Performance": memory :math:`32(J+1)GF` bytes, constant in :math:`C`.
  It has no :math:`S`, so its time is about :math:`(A + B)/C`, and it keeps
  gaining from more cores.

At most one core per geometry shares :math:`A`, and for "High Performance" at
most one core per element shares :math:`B`.  With a single geometry in all of
*measurement_sensors*, no kernel shares :math:`A` among the cores, though the
methods of the agenda may use them themselves, and the time is about
:math:`A + B` for "Low Memory" and "Indexed" and :math:`A + B/C` for "High
Performance".  The working memory of the agenda grows by one agenda call per
core for all kernels.  Memory bandwidth, which the cores share, is left out.

For example, take an antenna that samples 381 geometries per placement, at 13
placements, with 20 channels of 10 frequencies, so :math:`G = 4953`,
:math:`F = 200`, :math:`E = 260`, and :math:`P = 260 \cdot 381 = 99060`.  With
:math:`J = 500` and :math:`C = 16`:

- *measurement_jac*: :math:`8 \cdot 500 \cdot 260` bytes, about 1 MB.
- "Low Memory": :math:`16 \cdot 32 \cdot 501 \cdot 200` bytes, about 51 MB.
- "High Performance": :math:`32 \cdot 501 \cdot 4953 \cdot 200` bytes, about 15.9 GB.
- "Indexed": the 51 MB of "Low Memory", about 2.4 MB of list, and
  :math:`16 \cdot 8 \cdot 501 \cdot 20` bytes, about 1.3 MB, of contributions,
  as 20 elements have weights at each geometry.

Each further core adds about 3.2 MB to "Low Memory", 3.3 MB to "Indexed", and
nothing to "High Performance".
)",
      .values_and_desc = {Value{"LowMem", "Low Memory", "Use the low memory kernel"},
                          Value{"HighPerf", "High Performance", "Use the high performance kernel"},
                          Value{"Indexed", "Indexed", "Use the indexed kernel"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "InterpolationExtrapolation",
      .desc =
          R"(Instructions about how to handle extrapolation of interpolated data.
)",
      .values_and_desc = {Value{"None", "Do not allow interpolation of data outside the grid"},
                          Value{"Zero", "Values outside the grid are set to zero"},
                          Value{"Nearest", "Values outside the grid are set to the nearest grid point"},
                          Value{"Linear", "Linearly extrapolate the data outside the grid"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "PartitionFunctionsType",
      .desc =
          R"(Type of partition function data.

This is only used by the partition functions setup during compilation.
Each type represent a different way to store and compute the partition function of
an isotopologue in the :class:`~pyarts3.arts.PartitionFunctionsData` type.

The types, how they compute, and how they store date are as follows:

.. rubric:: Interp

The data contains the temperature grid and the corresponding partition function values.  The first
column is the temperature grid.  It must be in increasing order.
The second column is the partition function values.  The partition function is simply computed
by linear interpolation between the grid points.

.. rubric:: Constant

The data contains a single value, which is the constant partition function value for all temperatures.

.. rubric:: Coeff

The data contains the coefficients of a polynomial expansion
of the partition function as a function of temperature.
The partition function is computed by evaluating the polynomial at the given temperature.
The order of the polynomial is that the first coefficient is the constant term,
the second coefficient is the linear term, etc.

.. rubric:: StaticInterp

This is very similar to the interp-type, but the grid-values are assumed to be separated
by a constant value.  This constant value allows for a more efficient lookup and interpolation
of the partition function values.  The layout of the data is otherwise the same as the interp-type,
with the first column being the temperature grid and the second column being
the partition function values.  The first column must be in increasing order,
and the difference between consecutive values must be constant.

.. rubric:: LogInterp

The layout of the data is the same as for the interp-type, but the logarithm of the
partition function is interpolated linearly in the logarithm of the temperature.
Between two grid points, the partition function is thus a power law of the temperature.
This reproduces the rotational partition function of a rigid rotor, which goes as
:math:`T` or :math:`T^{3/2}`, also on a coarse grid, such as the one of the JPL catalog.
Outside the grid, the power law of the nearest pair of grid points is used.
All temperatures and partition function values must be positive.
)",
      .values_and_desc = {Value{"Interp", "Interpolate the data"},
                          Value{"Constant", "Use constant value"},
                          Value{"Coeff", "Use as polynomial coefficients"},
                          Value{"StaticInterp", "Interpolate the data, the temperature grid is known at compile time."},
                          Value{"LogInterp", "Interpolate the logarithm of the data in the logarithm of temperature"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "SpeciesEnum",
      .desc =
          R"(The valid species for the ARTS system.

This identifies the species that can be used in the ARTS system.
It has a wide variety of applications internally.

In terms of atmospheric properties, it is used to identify the VMR
in the atmosphere of the species.  See *AtmField* and *AtmPoint* for its usage there.

For some applications, it is simply an identifier in a map, such as for *SpeciesEnumVectors*.

For some applications, it is an extension of the species identifier, such
as for *QuantumIdentifier*, *QuantumLevelIdentifier*, and *SpeciesTag*.

The type is also used extensively in the absorption system for identifying species.
This includes *AbsorptionLookupTables*, *CIARecord*, *XsecRecord*, and *AbsorptionBands*.
Specifically for the latter, the inner *AbsorptionLine* type has line shape
parameters that are mapped to the species identifier.
)",
      .values_and_desc = {Value{"Bath", "AIR", "An all encompassing species for the rest of the atmosphere"}},
      .preferred_print = 1,
  });

  const auto& species = get_species_enum_info();
  for (const auto& info : species) {
    opts.back().values_and_desc.emplace_back(
        Value{info.longname, info.shortname, arts::format("A species short-named {}", info.shortname)});
  }
  opts.back().values_and_desc.emplace_back(
      Value{"unused",
            "unused",
            "This indicates that the concept of a species is not used in the context of the relevant option"});

  opts.emplace_back(EnumeratedOption{
      .name            = "QuantumNumberType",
      .desc            = R"(The type of value for a quantum number.
)",
      .values_and_desc = {Value{"alpha", "Quantum number \"alpha\" - Not in VAMDC"},
                          Value{"config", "Quantum number \"config\" - Not in VAMDC"},
                          Value{"ElecStateLabel", "Quantum number \"ElecStateLabel\""},
                          Value{"F", "Quantum number \"F\""},
                          Value{"F1", "Quantum number \"F1\""},
                          Value{"F10", "Quantum number \"F10\""},
                          Value{"F11", "Quantum number \"F11\""},
                          Value{"F12", "Quantum number \"F12\""},
                          Value{"F2", "Quantum number \"F2\""},
                          Value{"F3", "Quantum number \"F3\""},
                          Value{"F4", "Quantum number \"F4\""},
                          Value{"F5", "Quantum number \"F5\""},
                          Value{"F6", "Quantum number \"F6\""},
                          Value{"F7", "Quantum number \"F7\""},
                          Value{"F8", "Quantum number \"F8\""},
                          Value{"F9", "Quantum number \"F9\""},
                          Value{"I", "Quantum number \"I\""},
                          Value{"J", "Quantum number \"J\""},
                          Value{"K", "Quantum number \"K\""},
                          Value{"Ka", "Quantum number \"Ka\""},
                          Value{"Kc", "Quantum number \"Kc\""},
                          Value{"L", "Quantum number \"L\" - Not in VAMDC"},
                          Value{"Lambda", "Quantum number \"Lambda\""},
                          Value{"N", "Quantum number \"N\""},
                          Value{"Omega", "Quantum number \"Omega\""},
                          Value{"S", "Quantum number \"S\""},
                          Value{"Sigma", "Quantum number \"Sigma\""},
                          Value{"SpinComponentLabel", "Quantum number \"SpinComponentLabel\""},
                          Value{"asSym", "Quantum number \"asSym\""},
                          Value{"elecInv", "Quantum number \"elecInv\""},
                          Value{"elecRefl", "Quantum number \"elecRefl\""},
                          Value{"elecSym", "Quantum number \"elecSym\""},
                          Value{"kronigParity", "Quantum number \"kronigParity\""},
                          Value{"l", "Quantum number \"l\""},
                          Value{"l1", "Quantum number \"l1\""},
                          Value{"l10", "Quantum number \"l10\""},
                          Value{"l11", "Quantum number \"l11\""},
                          Value{"l12", "Quantum number \"l12\""},
                          Value{"l2", "Quantum number \"l2\""},
                          Value{"l3", "Quantum number \"l3\""},
                          Value{"l4", "Quantum number \"l4\""},
                          Value{"l5", "Quantum number \"l5\""},
                          Value{"l6", "Quantum number \"l6\""},
                          Value{"l7", "Quantum number \"l7\""},
                          Value{"l8", "Quantum number \"l8\""},
                          Value{"l9", "Quantum number \"l9\""},
                          Value{"n", "Quantum number \"n\" - Not in VAMDC"},
                          Value{"parity", "Quantum number \"parity\""},
                          Value{"r", "Quantum number \"r\""},
                          Value{"rotSym", "Quantum number \"rotSym\""},
                          Value{"rovibSym", "Quantum number \"rovibSym\""},
                          Value{"sym", "Quantum number \"sym\""},
                          Value{"tau", "Quantum number \"tau\" - Not in VAMDC"},
                          Value{"term", "Quantum number \"term\" - Not in VAMDC"},
                          Value{"v", "Quantum number \"v\""},
                          Value{"v1", "Quantum number \"v1\""},
                          Value{"v10", "Quantum number \"v10\""},
                          Value{"v11", "Quantum number \"v11\""},
                          Value{"v12", "Quantum number \"v12\""},
                          Value{"v2", "Quantum number \"v2\""},
                          Value{"v3", "Quantum number \"v3\""},
                          Value{"v4", "Quantum number \"v4\""},
                          Value{"v5", "Quantum number \"v5\""},
                          Value{"v6", "Quantum number \"v6\""},
                          Value{"v7", "Quantum number \"v7\""},
                          Value{"v8", "Quantum number \"v8\""},
                          Value{"v9", "Quantum number \"v9\""},
                          Value{"vibInv", "Quantum number \"vibInv\""},
                          Value{"vibRefl", "Quantum number \"vibRefl\""},
                          Value{"vibSym", "Quantum number \"vibSym\""},
                          Value{"G", "Quantum number \"G\""},
                          Value{"U", "Quantum number \"U\""}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "IsoRatioOption",
      .desc            = R"(The type of isotopologue ratio to use.
)",
      .values_and_desc = {Value{"Builtin", "Use the built-in isotopologue ratio"},
                          Value{"Hitran", "Use the HITRAN isotopologue ratio"},
                          Value{"None", "Do not use an isotopologue ratio"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "LineShapeModelType",
      .desc = R"(The type of line shape model to use.
)",
      .values_and_desc =
          {Value{"T0", ":math:`X_0`"},
           Value{"T1", R"(:math:`X_0 \left(\frac{T_0}{T}\right) ^ {X_1}`)"},
           Value{
               "T2",
               R"(:math:`X_0 \left(\frac{T_0}{T}\right) ^ {X_1} \left[1 + X_2 \log\left(\frac{T_0}{T}\right)\right]`)"},
           Value{"T3", R"(:math:`X_0 + X_1 \left(T - T_0\right)`)"},
           Value{"T4",
                 R"(:math:`\left[X_0 + X_1 \left(\frac{T_0}{T} - 1\right)\right] \left(\frac{T_0}{T}\right)^{X_2}`)"},
           Value{"T5", R"(:math:`X_0 \left(\frac{T_0}{T}\right)^{\frac{1}{4} + \frac{3}{2}X_1}`)"},
           Value{
               "AER",
               R"(:math:`X(200) = X_0`; :math:`X(250) = X_1`; :math:`X(298) = X_2`; :math:`X(340) = X_3`;  Linear interpolation in between)"},
           Value{"DPL", R"(:math:`X_0 \left(\frac{T_0}{T}\right) ^ {X_1} + X_2 \left(\frac{T_0}{T}\right) ^ {X_3}`)"},
           Value{"POLY", R"(:math:`X_0 + X_1 T + X_2 T ^ 2 + X_3 T ^ 3 + \cdots`)"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "LineShapeModelCoefficient",
      .desc =
          R"(The type of line shape model coefficients.  See :class:`~pyarts3.arts.LineShapeModelType` for more information.
)",
      .values_and_desc =
          {Value{"X0", "x0", "X_0", "x_0", ":math:`X_0`"},
           Value{"X1", "x1", "X_1", "x_1", ":math:`X_1`"},
           Value{"X2", "x2", "X_2", "x_2", ":math:`X_2`"},
           Value{"X3", "x3", "X_3", "x_3", ":math:`X_3`"},

           Value{
               "unused", "invalid", "not-an-option", "void", "Unused placeholder - for default values of inactivity"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "LineShapeModelVariable",
      .desc            = R"(The type of line shape model variable.
)",
      .values_and_desc = {Value{"G0", "Pressure broadening speed-independent"},
                          Value{"D0", "Pressure f-shifting speed-dependent"},
                          Value{"G2", "Pressure broadening speed-dependent"},
                          Value{"D2", "Pressure f-shifting speed-independent"},
                          Value{"FVC", "Frequency of velocity-changing collisions"},
                          Value{"ETA", "Correlation"},
                          Value{"Y", "First order line mixing coefficient"},
                          Value{"G", "Second order line mixing coefficient"},
                          Value{"DV", "Second order line mixing f-shifting"},

                          Value{"unused", "Unused placeholder - for default values of inactivity"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "PathPositionType",
      .desc            = R"(A type of position in a path.
)",
      .values_and_desc = {Value{"atm", "Atmospheric position"},
                          Value{"space", "Space position"},
                          Value{"surface", "Surface position"},
                          Value{"unknown", "Unknown position"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "LineByLineVariable",
      .desc            = R"(A type of line by line variable.
)",
      .values_and_desc = {Value{"f0", "Central frequency [Hz]"},
                          Value{"e0", "Lower level energy [J]"},
                          Value{"a", "Einstein coefficient [1/s]"},

                          Value{"unused", "Unused placeholder - for default values of inactivity"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "LineByLineCutoffType",
      .desc            = R"(A type of line by line cutoff.
)",
      .values_and_desc = {Value{"None", "No cutoff"},
                          Value{"ByLine", "Line's are cut 1-by-1 around f0 in *AbsorptionLine*"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "LineByLineLineshape",
      .desc = R"(A type of line shape for line by line calculations.
)",
      .values_and_desc =
          {Value{"VP_LTE", "Voigt in local thermodynamic equilibrium"},
           Value{"VP_LTE_MIRROR", "Voigt in local thermodynamic equilibrium with negative frequency lines"},
           Value{"VP_LINE_NLTE", "Voigt in non-local thermodynamic equilibrium with level-by-level data"},
           Value{"VP_ECS_MAKAROV", "Voigt using Makarov's method of energy-corrected sudden for line mixing of O2"},
           Value{"VP_ECS_HARTMANN", "Voigt using Hartmann's method of energy-corrected sudden for line mixing of CO2"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "PlanetOrMoonType",
      .desc            = R"(The type of planetary body that should be considered.
)",
      .values_and_desc = {Value{"Earth", "Planet is Earth"},
                          Value{"Io", "Moon is Jupiter's Io"},
                          Value{"Jupiter", "Planet is Jupiter"},
                          Value{"Mars", "Planet is Mars"},
                          Value{"Venus", "Planet is Venus"},
                          Value{"Moon", "Moon is Earth's Moon"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "MissingFieldComponentError",
      .desc            = R"(What kind of error is it to miss a field component?
)",
      .values_and_desc = {Value{"Throw", "It is required to have the field component."},
                          Value{"Zero", "It is fine to zero-out the field component if it is missing."},
                          Value{"Ignore", "It is fine to ignore that the field component is missing."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "HydrostaticPressureOption",
      .desc =
          R"(These options control how the hydrostatic pressure is computed in ARTS.

There are two main options for how the hydrostatic pressure is computed:

- ``HydrostaticEquation``:
  
  .. math::
    p(z + h) = p(z) \left(1 - \frac{Mg}{RT} h\right) 

- ``HypsometricEquation``:
  
  .. math::
    p(z + h) = p(z) e^{- \frac{Mg}{RT} h } 

where:

- :math:`p(z)` is the pressure at altitude :math:`z`,
- :math:`h` is the altitude difference,
- :math:`g` is the gravitational acceleration at altitude :math:`z`,
- :math:`M` is the mean mass of an atmospheric molecule at altitude :math:`z`,
- :math:`R` is the specific gas constant, and
- :math:`T` is the temperature at altitude :math:`z`.

.. note::
  We do not consider The Eötvös effect.
  We also do not directly support "virtual temperature", as water is baked into :math:`M`.
  Either addition would be welcome as additional options.
)",
      .values_and_desc = {Value{"HydrostaticEquation", "Piece-wise linearly decaying with altitude."},
                          Value{"HypsometricEquation", "Piece-wise exponentially decaying with altitude."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "FieldComponent",
      .desc            = R"(Selection of a field component
)",
      .values_and_desc = {Value{"u", "U", "East component"},
                          Value{"v", "V", "North component"},
                          Value{"w", "W", "Up component"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "SpectralRadianceUnitType",
      .desc            = R"(Choice of spectral radiance unit in conversions.

The conversion to one of these units is expected to happen from the internal representation of *spectral_rad*,
which is [W / m :math:`^{2}` Hz sr].

For the description below, :math:`c` is the speed of light, :math:`k` is the Boltzmann constant,
:math:`f` is the frequency, :math:`h` is the Planck constant, :math:`F(x)` is the conversion
function, and :math:`[I,\; Q,\; U,\; V]` is the internal representation of *spectral_rad*.

For Rayleigh-Jeans brightness temperature the conversion is:

.. math::

    [I,\; Q,\; U,\; V] \rightarrow [F(I),\; F(Q),\; F(U),\; F(V)]

.. math::

    F(x) = \frac{c^2 x}{2 k f^2}.

For Planck brightness temperature the conversion is:

.. math::

    [I,\; Q,\; U,\; V] \rightarrow
        \left[
        F(I),\;
        F\left( \frac{I + Q}{2} \right) - F\left( \frac{I - Q}{2} \right),\;
        F\left( \frac{I + U}{2} \right) - F\left( \frac{I - U}{2} \right),\;
        F\left( \frac{I + V}{2} \right) - F\left( \frac{I - V}{2} \right)
        \right]

.. math::

    F(x) = \frac{h f}{k}\left[\log\left(1 + \frac{2hf^3}{c^2x}\right)\right]^{-1}.

The spectral radiance per wavelength [W / m :math:`^{2}` m sr] is:

.. math::

    [I,\; Q,\; U,\; V] \rightarrow [F(I),\; F(Q),\; F(U),\; F(V)]

.. math::

    F(x) = \frac{f^2 x}{c}.

The spectral radiance per wavenumber [W / m :math:`^{2}` m :math:`^{-1}` sr] is:

.. math::

    [I,\; Q,\; U,\; V] \rightarrow [F(I),\; F(Q),\; F(U),\; F(V)]

.. math::

    F(x) = \frac{x}{c}.

Lastly, the unit option of course just retains the current state [W / m :math:`^{2}` Hz sr].

.. math::

    [I,\; Q,\; U,\; V] \rightarrow [I,\; Q,\; U,\; V].
)",
      .values_and_desc = {Value{"RJBT",
                                "Tr",
                                R"(Rayleigh-Jeans brightness temperature [K]
)"},
                          Value{"PlanckBT",
                                "Tb",
                                R"(Planck brightness temperature [K]
)"},
                          Value{"W_m2_m_sr", "W/(m^2 m sr)", "Spectral radiance wavelength [W / m :math:`^{2}` m sr]"},
                          Value{"W_m2_m1_sr",
                                "W/(m^2 m-1 sr)",
                                "Spectral radiance wavenumber [W / m :math:`^{2}` m :math:`^{-1}` sr]"},
                          Value{"unit", "1", "Unit spectral radiance [W / m :math:`^{2}` Hz sr]"}},
  });

  //  Manually defined to decouple file IO from this code (to help species generation)
  //   opts.emplace_back(EnumeratedOption{
  //       .name = "FileType",
  //       .desc = "A choice of file format types.\n",
  //       .values_and_desc =
  //           {Value{"ascii", "ASCII", "Ascii", "text", "Save as ASCII"},
  //            Value{"zascii", "ZASCII", "Zip", "zip", "Save as zipped ASCII"},
  //            Value{"binary", "BINARY", "Binary", "bin", "Save as binary data"}},
  //   });

  opts.emplace_back(EnumeratedOption{
      .name            = "EarthEllipsoid",
      .desc            = "Choice of ellipsoid.\n",
      .values_and_desc = {Value{"Sphere", R"(A spherical Earth.)"},
                          Value{"WGS84", R"(The reference ellipsoid used by the GPS system.)"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "IoEllipsoid",
      .desc            = "Choice of ellipsoid.\n",
      .values_and_desc = {Value{
          "Sphere", "A spherical planetesimal.  The radius is taken from a report of the IAU/IAG Working Group."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "EuropaEllipsoid",
      .desc            = "Choice of ellipsoid.\n",
      .values_and_desc = {Value{"Sphere", "A spherical planetesimal."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "GanymedeEllipsoid",
      .desc            = "Choice of ellipsoid.\n",
      .values_and_desc = {Value{"Sphere", "A spherical planetesimal."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "JupiterEllipsoid",
      .desc = "Choice of ellipsoid.\n",
      .values_and_desc =
          {Value{"Sphere", "A spherical planet. The radius is taken from a report of the IAU/IAG Working Group."},
           Value{"Ellipsoid",
                 "A reference ellipsoid with parameters taken from a report of the IAU/IAG Working Group."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "MarsEllipsoid",
      .desc = "Choice of ellipsoid.\n",
      .values_and_desc =
          {Value{"Sphere", "A spherical planet. The radius is taken from a report of the IAU/IAG Working Group."},
           Value{"Ellipsoid",
                 "A reference ellipsoid with parameters taken from a report of the IAU/IAG Working Group."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "MoonEllipsoid",
      .desc = "Choice of ellipsoid.\n",
      .values_and_desc =
          {Value{"Sphere", "A spherical planet. The radius is taken from a report of the IAU/IAG Working Group."},
           Value{"Ellipsoid",
                 "A reference ellipsoid with parameters taken from a report of the IAU/IAG Working Group."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "VenusEllipsoid",
      .desc            = "Choice of ellipsoid.\n",
      .values_and_desc = {Value{"Sphere",
                                "A spherical planet. The radius is taken from a report of the IAU/IAG Working Group."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "SpeciesTagType",
      .desc            = "Type of species tag.\n",
      .values_and_desc = {Value{"Plain", "A plain species, one isotopologue or all."},
                          Value{"Predefined", "A predefined model species."},
                          Value{"Cia", "A pair of collision-induced species."},
                          Value{"XsecFit", "A cross-section fitting of a model species."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "PolarizationChoice",
      .desc = R"(Named polarization states to help create relevant *Stokvec*.

.. note::
    These are just user friendly suggestions and it is recommended to
    create the correct *Stokvec* manually if the desired polarization state is not
    represented below.

    Also, be aware that the unit of, e.g., *spectral_rad* (often the last choice
    of *SpectralRadianceUnitType*) is in Kelvin, then the code below will give the
    measured brightness temperature in Kelvin for these polarization states, but that
    if the unit is still in Watts of any kind, then the code below will give 2 times
    the select polarized brightness temperatures (but still the correct unpolarized
    radiation).
)",
      .values_and_desc =
          {Value{"I", "I", "No polarization state [1, 0, 0, 0]."},
           Value{"Q", "Q", "Difference between vertical and horizontal linear polarization [0, 1, 0, 0]."},
           Value{"U", "U", "Difference between plus and minus 45 degrees linear polarization [0, 0, 1, 0]."},
           Value{"V", "V", "Difference between right and left circular polarization [0, 0, 0, 1]."},
           Value{"Iv", "Iv", "Vertical linear polarization [1, 1, 0, 0]."},
           Value{"Ih", "Ih", "Horizontal linear polarization [1, -1, 0, 0]."},
           Value{"Ip45", "I+45", "Plus 45 degrees linear polarization [1, 0, 1, 0]."},
           Value{"Im45", "I-45", "Minus 45 degrees linear polarization [1, 0, -1, 0]."},
           Value{"Ilhc", "LC", "Left circular polarization [1, 0, 0, -1]."},
           Value{"Irhc", "RC", "Right circular polarization [1, 0, 0, 1]."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "TMatrixSizeDistribution",
      .desc =
          R"(The size distributions of the T-matrix solver's random-orientation average.

They are those of the original solver, with the parameters a, b and gamma.
See :doc:`concept.tmatrix` for the definitions.
)",
      .values_and_desc = {Value{"modified_gamma", "The modified gamma distribution"},
                          Value{"lognormal", "The log-normal distribution"},
                          Value{"power_law", "The power law of Hansen and Travis (1974)"},
                          Value{"gamma", "The gamma distribution"},
                          Value{"modified_power_law", "The modified power law"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "ParticleShape",
      .desc =
          R"(Rotationally symmetric particle shapes, as used by the T-matrix solver.

An aspect ratio, and for Chebyshev particles a degree, complete the shape.
See :doc:`concept.tmatrix` for the definitions.
)",
      .values_and_desc =
          {Value{"spheroid", "A spheroid; the aspect ratio is its horizontal-to-rotational axis ratio"},
           Value{"cylinder", "A finite circular cylinder; the aspect ratio is its diameter-to-length ratio"},
           Value{
               "chebyshev",
               R"(A Chebyshev particle :math:`r_0 (1 + e \cos(n \theta))` of degree :math:`n \geq 1`; the aspect ratio is its deformation :math:`e`, :math:`|e| < 1`)"},
           Value{
               "generalized_chebyshev",
               R"(A generalized Chebyshev particle :math:`r_0 (1 + c_0 + \sum_n c_n \cos(n \theta))`, e.g. a distorted water drop; its coefficients :math:`c_0, \ldots, c_N` are given, and the aspect ratio is not used)"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "EquivalentRadius",
      .desc =
          R"(The sphere whose radius gives the size of a nonspherical particle, as used by the T-matrix solver.

See :doc:`concept.tmatrix` for the conversion.
)",
      .values_and_desc = {Value{"volume", "The radius of the sphere of the same volume"},
                          Value{"surface_area", "The radius of the sphere of the same surface area"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "ParticulateProperty",
      .desc =
          R"(Numerical properties used to numerically represent particle populations.
)",
      .values_and_desc = {Value{"MassDensity", "m", "Mass density in kg/m^{-3}"},
                          Value{"NumberDensity", "n", "Number density in m^{-3}"},
                          Value{"DMax", "dmax", "Maximum particle diameter in m."},
                          Value{"DVeq", "dveq", "Volume-equivalent diameter in m."},
                          Value{"Extinction", "ext", "Extinction in m^{-1}"},
                          Value{"SingleScatteringAlbedo", "ssa", "Single scattering albedo"},
                          Value{"ShapeParameter", "ShapeParameter", "PSD shape parameter in arbitrary units."},
                          Value{"InterceptParameter",
                                "InterceptParameter",
                                "PSD intercept parameter in arbitrary units."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "SizeParameter",
      .desc =
          R"(Parameters used to represent the size of particles.
)",
      .values_and_desc = {Value{"Mass", "m", "Mass in kg"},
                          Value{"DMax", "d_max", "Maximum diameter in m"},
                          Value{"DVeq", "d_veq", "Volume-equivalent diameter in m"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "TransmittanceOption",
      .desc = R"(A flag for the transmittance model for a layer.
      
.. rubric:: Constant 

Classical and simple solution for radiative transfer over a layer with constant
absorption coefficient and constant source function.  Strong variations in the
source function over the layer may lead to inaccuracies.  As may strong variations
in the absorption coefficient.  Easiest and fastest method.

.. math::

    K = \frac{K_0 + K_1}{2}

    T = e^{- K r}

    J = \frac{J_0 + J_1}{2}

    I_1 = J + T (I_0 - J)

.. rubric:: Linear in Tau 
      
The source is linearly changing throughout the layer, while the absorption coefficient
is constant.  Has better accuracy than the ``constant`` option when the source function
varies strongly over the layer.  Still fast to compute, but the linear evolution operator
is highly complex, meaning we do not know if there are edge cases where this method fails.

.. math::

    K = \frac{K_0 + K_1}{2}

    T = e^{- K r}

    J = J_0 + \left(J_1 - J_0\right) \frac{s}{r}

    \Lambda_0 = \frac{1}{r}K^{-1} \left(1 - T\right)

    I_1 = J_1 + T (I_0 - J_0) + \Lambda_0 \left(J_0 - J_1\right)

.. rubric:: Linear in Propagation

The source is linearly changing throughout the layer, as is the the absorption coefficient.
Has theoretically better accuracy than the ``constant`` or ``linsrc`` options when the layer
varies strongly.  However, the linear evolution operator is highly complex, meaning we do not know
if there are edge cases where this method fails, and comparisons with ``linsrc`` have not shown
significant improvements in accuracy - in fact sometimes it is less accurate.

.. math::

    K = K_0 + \left(K_1 - K_0\right) \frac{s}{r}

    T = e^{- K r}

    J = J_0 + \left(J_1 - J_0\right) \frac{s}{r}

    \Lambda_0 = \frac{1}{r}\sqrt{\frac{K_1 - K_0}{2r}} ^ {-1} \left[D\left(\frac{1}{2} \sqrt{\frac{K_1 - K_0}{2r}} ^ {-1} K_1\right) - T D\left(\frac{1}{2} \sqrt{\frac{K_1 - K_0}{2r}} ^ {-1} K_0\right)\right],\;\mathrm{when}\;K_1 \gt K_0

    \Lambda_0 = \frac{1}{r}K^{-1} \left(1 - T\right),\; \mathrm{otherwise}

    I_1 = J_1 + T (I_0 - J_0) + \Lambda_0 \left(J_0 - J_1\right)

where :math:`D(x)` is the Dawson function.

.. rubric:: Magnus Operator

Ordered expansion of the matrix exponential of the linear evolution operator.
Has theoretically better accuracy than the ``constant`` when the layer varies strongly.
This is experimental code and might not work in all cases, especially as :math:`r` grows
large, considering the quadratic term.
Here :math:`K_0` is at the start of the propagation step and :math:`K_1` is at its end.

.. math::

    O = r\frac{K_0 + K_1}{2} - \frac{r^2}{12} \left[K_1, K_0\right]

    T = e^{- O}

    J = \frac{J_0 + J_1}{2}

    I_1 = J + T (I_0 - J)

For a linearly varying source, the augmented Magnus system gives

.. math::

    L = \phi_1(-O)\left[1 - \frac{r}{12}(K_1-K_0)\right]

    \phi_1(X) = X^{-1}\left(e^X-1\right)

    I_1 = J_1 + T(I_0-J_0) + L(J_0-J_1)
)",
      .values_and_desc =
          {Value{"constant", "linear", "lin", R"(Linear transmittance step over constant layer)"},
           Value{"linsrc",
                 "lintau",
                 "linear-in-tau",
                 R"(Linear transmittance step over constant absorption layer with linear source function)"},
           Value{"linprop",
                 "linsrcprop",
                 "linear-in-propagation",
                 R"(Linear transmittance step over non-constant absorption layer with linear source function)"},
           Value{
               "magop",
               "MagnusOperator",
               "magnus",
               R"(Magnus operator transmittance step over a linearly varying propagation matrix with layer-average source)"},
           Value{"magop_linsrc",
                 "MagnusLinearSource",
                 "magnus-linear-source",
                 R"(Magnus operator transmittance step over a linearly varying propagation matrix and source)"}},
  });

  opts.emplace_back(EnumeratedOption{
      .name = "WorkspaceInitialization",
      .desc = "A flag for how initialize a new workspace.\n",
      .values_and_desc =
          {Value{"FromGlobalDefaults",
                 "Initialize from the global default workspace - these are the variables that have defaults."},
           Value{"Empty", "Initialize an empty workspace - no variables are set."}},
  });

  opts.emplace_back(EnumeratedOption{
      .name            = "ZeemanPolarization",
      .desc            = "A flag for the Zeeman polarization state.\n",
      .values_and_desc = {Value{"sm", "Parallel minus."},
                          Value{"pi", "Perpendicular."},
                          Value{"sp", "Parallel plus."},
                          Value{"no", "Unpolarized."}},
  });

  fix_static(opts);
  return opts;
}
}  // namespace

const std::vector<EnumeratedOption>& internal_options() {
  static std::vector<EnumeratedOption> opts = internal_options_create();
  return opts;
}

std::string_view EnumeratedOption::sz() const {
  const auto n = values_and_desc.size();

  if (n <= 2) return "bool";
  if (n <= std::numeric_limits<unsigned char>::max()) return "unsigned char";
  if (n <= std::numeric_limits<unsigned short>::max()) return "unsigned short";
  if (n <= std::numeric_limits<unsigned int>::max()) return "unsigned int";
  if (n <= std::numeric_limits<unsigned long>::max()) return "unsigned long";

  throw std::runtime_error("Too many values for enum class");
}

namespace {
std::string add_spaces(const std::string& s, int numspaces) {
  std::string       out{s};
  const std::string newline{"\n"};
  const std::string spaces = std::string(numspaces, ' ') + '\n';
  replace(out, newline, spaces);
  return out;
}
}  // namespace

std::string EnumeratedOption::docs() const {
  std::string out;

  if (values_and_desc.empty()) { throw std::runtime_error("No values for enum class " + name); }

  const auto n = values_and_desc.front().size();

  out += arts::format("{}\n\nValid options:\n\n", desc);
  for (auto& v : values_and_desc) {
    if (v.empty()) { throw std::runtime_error("Empty value for enum class " + name); }

    std::string_view x = "- "sv;
    for (auto& s : v | stdv::take(n - 1)) { out += arts::format(R"({}``"{}"``)", std::exchange(x, " or "sv), s); }
    out += arts::format(":\n  {}\n", add_spaces(v.back(), 2));
  }

  return out;
}

std::string EnumeratedOption::tail() const {
  std::string out;

  if (values_and_desc.empty()) { throw std::runtime_error("No values for enum class " + name); }

  const auto m = values_and_desc.size();
  const auto n = values_and_desc.front().size();

  // Create good enum function
  out += arts::format(
      R"(
template<> constexpr bool good_enum<{0}>({0} x) noexcept {{
    const auto v = static_cast<std::size_t>(x);
    return v < {1};
}}

template<> struct enumdocs<{0}> {{
  static std::string_view str() noexcept;
  static constexpr std::string_view name = "{0}"sv;
}};

namespace enumtyps {{
    inline constexpr std::array {0}Types = {{)",
      name,
      m);

  for (const auto& v : values_and_desc) out += arts::format("{}::{}, ", name, v[0]);

  out += R"(};
}  // namespace enumtyps

namespace enumstrs {
)";

  // Create enumstrs arrays
  for (std::size_t i = 0; i < n - 1; i++) {
    out += arts::format(
        R"(
  template <> struct enum_str_data<{}, {}> {{
    static constexpr std::array strs={{)",
        name,
        i);

    for (const auto& v : values_and_desc) out += arts::format(R"("{}"sv, )", v[i]);

    out += "};\n  };\n";
  }

  out += arts::format(
      R"(
  template <int i={0}>
  inline constexpr auto {1}Names = enum_str_data<{1}, i>::strs;
}}  // namespace enumstrs

template <int i={0}> constexpr std::string_view toString({1} x)  requires(i >= 0 and i < {2}) {{
  if (good_enum(x))
    return enumstrs::{1}Names<i>[static_cast<std::size_t>(x)];
  return "BAD {1}"sv;
}}

template<> constexpr {1} to<{1}>(const std::string_view x) {{
  switch (enumhash(x)) {{)",
      preferred_print,
      name,
      n - 1);

  /* Emit one case per distinct key hash.  head() guarantees that a key
   * maps to at most one enumerator; keys are visited in name-array
   * order and inserted once, preserving the old first-match-wins scan
   * order.  Keys sharing a hash are chained under a single case label;
   * the full-string comparison decides, so collisions are harmless. */
  {
    std::vector<std::pair<std::string_view, std::string_view>> keys;
    std::set<std::string_view>                                 seen;
    for (std::size_t i = 0; i < n - 1; i++) {
      for (const auto& v : values_and_desc) {
        if (seen.insert(v[i]).second) keys.emplace_back(v[i], v[0]);
      }
    }

    std::vector<std::vector<std::size_t>> groups;
    for (std::size_t k = 0; k < keys.size(); k++) {
      const auto h = enumhash(keys[k].first);
      auto       g = stdr::find_if(groups, [&](const auto& group) { return enumhash(keys[group.front()].first) == h; });
      if (g == groups.end()) {
        groups.push_back({k});
      } else {
        g->push_back(k);
      }
    }

    for (const auto& group : groups) {
      out += arts::format(R"(
    case enumhash("{0}"sv):)",
                          keys[group.front()].first);
      for (const auto k : group) {
        out += arts::format(R"(
      if (x == "{0}"sv) return {1}::{2};)",
                            keys[k].first,
                            name,
                            keys[k].second);
      }
      out += "\n      break;";
    }
  }

  out += "\n  }\n";
  out += arts::format(
      R"(
  throw std::runtime_error(arts::format(R"-x-(Bad input "{{}}"

See https://atmtools.github.io/arts-docs-master/pyarts3.arts.{0}.html for valid options.
)-x-", x));
}}

namespace enumsize {{ inline constexpr std::size_t {0}Size = {2}; }}

template<> struct std::formatter<{0}> {{
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() {{ return *this; }}
  [[nodiscard]] constexpr auto& inner_fmt() const {{ return *this; }}

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {{
    return parse_format_tags(tags, ctx);
  }}

  template <class FmtContext>
  FmtContext::iterator format(const {0}& v, FmtContext& ctx) const {{
    return tags.format(ctx, tags.quote(), toString<{3}>(v), tags.quote());
  }}
}};)",
      name,
      docs(),
      m,
      preferred_print);

  return out;
}

std::string EnumeratedOption::head() const {
  static constexpr std::array reserved = {"get_options", "get_options_as_string"};

  for (auto& v : values_and_desc) {
    for (auto& s : v) {
      if (stdr::any_of(reserved, [&s](auto& r) { return r == s; })) {
        throw std::runtime_error("Reserved word used for enum class " + name);
      }
    }
  }

  if (values_and_desc.empty()) { throw std::runtime_error("No values for enum class " + name); }

  [*this]() {
    std::vector<std::string> all;
    all.reserve(values_and_desc.size() * (values_and_desc.front().size()));
    for (auto v : values_and_desc) {
      stdr::sort(v);
      auto last = std::unique(v.begin(), v.end());
      all.insert(all.end(), v.begin(), last);
    }
    stdr::sort(all);
    auto last = std::unique(all.begin(), all.end());
    if (last != all.end()) {
      throw std::runtime_error("Duplicate value for enum class \"" + name + "\" \"" + *last + "\"");
    }
  }();

  const auto n = values_and_desc.front().size();

  if (n < 2) {
    throw std::runtime_error(
        "No descriptions for enum class " + name +
        std::string{".\nThe Value should be {main_key, alt_keys..., description}, where alt_keys... can be omitted"});
  }

  std::string out;

  out += arts::format("enum class {} : {} {{\n", name, sz());

  for (const auto& v : values_and_desc) {
    if (v.empty()) { throw std::runtime_error("Empty value for enum class " + name); }

    if (v.size() != n) { throw std::runtime_error("Mismatched value count for enum class " + name + " " + v.front()); }

    out += arts::format("{}, ", v[0]);
  }

  out += "};\n\n";

  return out;
}

std::string EnumeratedOption::impl() const {
  std::string out;

  out += arts::format(
      R"(
std::string_view enumdocs<{0}>::str() noexcept {{
  return R"-ENUMDOC-({1})-ENUMDOC-"sv;
}}

static_assert(arts_formattable<{0}>, "Not good: {0}");
)",
      name,
      docs());

  return out;
}
