"""Radiance, flux and heating-rate diagnostics in SI units.

Positive heating warms the gas. Flux is positive upward. These functions
accept NumPy-compatible arrays, including ARTS arrays. They do not infer
brightness-temperature units, density, gravity or heat capacity.
"""

import numpy as np


def _finite(value, name):
    value = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must be finite")
    return value


def _positive(value, name):
    value = _finite(value, name)
    if np.any(value <= 0):
        raise ValueError(f"{name} must be positive")
    return value


def integrate_spectral(values, frequency, *, axis=0, weights=None):
    """Integrate per-Hz values over frequency [Hz].

    Use trapezoidal integration on an increasing frequency grid, or explicit
    quadrature weights [Hz]. A single frequency requires an explicit weight.
    The selected axis is removed; all other axes retain their order.
    """
    values = _finite(values, "values")
    frequency = _positive(frequency, "frequency")
    if frequency.ndim != 1 or frequency.size == 0 or np.any(np.diff(frequency) <= 0):
        raise ValueError("frequency must be a nonempty increasing vector")
    values = np.moveaxis(values, axis, 0)
    if values.shape[0] != frequency.size:
        raise ValueError("frequency does not match the spectral axis")
    if weights is None:
        if frequency.size < 2:
            raise ValueError("one frequency requires an explicit integration weight")
        return np.trapezoid(values, frequency, axis=0)
    weights = _finite(weights, "weights")
    if weights.shape != frequency.shape:
        raise ValueError("weights must match frequency")
    return np.einsum("f,f...->...", weights, values)


def flux_from_radiance(radiance, zenith, zenith_weights, *, axis=-1):
    """Return positive upward and downward hemispheric irradiance.

    ``radiance`` is Stokes I in W/(m² sr), or W/(m² sr Hz) for spectral
    data. It must be azimuth independent or already averaged over azimuth.
    ``zenith`` [degrees] uses the ARTS viewing direction: 0 looks upward
    (downward travelling radiation), 180 looks downward. ``zenith_weights``
    integrate d(cos(zenith)), e.g. double-Gauss weights summing to two.
    This is a projected flux integral, with a cosine factor, not actinic flux.
    The angular axis is removed. Select Stokes I before calling this function.
    """
    radiance = np.moveaxis(_finite(radiance, "radiance"), axis, -1)
    zenith = _finite(zenith, "zenith")
    weights = _finite(zenith_weights, "zenith_weights")
    if zenith.ndim != 1 or zenith.size == 0 or weights.shape != zenith.shape:
        raise ValueError("zenith and zenith_weights must be matching vectors")
    if np.any((zenith < 0) | (zenith > 180)) or np.any(weights < 0):
        raise ValueError("zenith must be in [0, 180] and weights nonnegative")
    if radiance.shape[-1] != zenith.size:
        raise ValueError("zenith does not match the angular axis")
    mu = np.cos(np.deg2rad(zenith))
    up = np.einsum("...a,a->...", radiance, 2 * np.pi * weights * np.maximum(-mu, 0))
    down = np.einsum("...a,a->...", radiance, 2 * np.pi * weights * np.maximum(mu, 0))
    return up, down


def from_flux_divergence(pressure_derivative, heat_capacity, gravity):
    """Convert d(F_up-F_down)/dp [W m⁻² Pa⁻¹] to heating [K/s].

    Hydrostatic balance gives H = g/cp * dF/dp. ``heat_capacity`` is mass
    specific [J/(kg K)] and ``gravity`` is positive [m/s²]. Scalar or array
    inputs follow NumPy broadcasting. No boundary values are suppressed.
    """
    return (
        _finite(pressure_derivative, "pressure_derivative")
        * _positive(gravity, "gravity")
        / _positive(heat_capacity, "heat_capacity")
    )


def from_flux(net_upward_flux, pressure, heat_capacity, gravity, *, axis=-1):
    """Derive heating [K/s] from upward net flux [W/m²] at pressure levels.

    ``pressure`` [Pa] is a strictly monotonic vector of at least three levels,
    in either order. Differentiate with the nonuniform three-point formula,
    including second-order one-sided boundaries. Heat capacity [J/(kg K)]
    and gravity [m/s²] broadcast against the result, whose shape is unchanged.
    For spectral flux, integrate over frequency first, or integrate the
    resulting spectral heating rates afterwards if g/cp is frequency independent.
    """
    flux = _finite(net_upward_flux, "net_upward_flux")
    pressure = _positive(pressure, "pressure")
    dp = np.diff(pressure) if pressure.ndim == 1 else np.array([])
    if (
        pressure.ndim != 1
        or pressure.size < 3
        or not (np.all(dp > 0) or np.all(dp < 0))
    ):
        raise ValueError(
            "pressure must be a strictly monotonic vector of at least three levels"
        )
    if flux.shape[axis] != pressure.size:
        raise ValueError("pressure does not match the profile axis")
    return from_flux_divergence(
        np.gradient(flux, pressure, axis=axis, edge_order=2), heat_capacity, gravity
    )


def from_optical_depth_derivative(dfdt, extinction, density, heat_capacity):
    """Convert DISORT/VDISORT DFDT to heating, preserving spectral axes.

    DFDT is d(F_up-F_down)/d(tau), with optical depth increasing downward.
    Multiply by physical extinction [1/m], then divide by density [kg/m³]
    and mass specific heat capacity [J/(kg K)]. The result is K/s for
    broadband inputs or K/(s Hz) for spectral inputs. All quantities must
    refer to the same locations; array inputs follow NumPy broadcasting.
    Extinction must be total extinction, not absorption (DFDT already
    contains the single-scattering-albedo factor). If a low-level VDISORT
    calculation uses externally rescaled optical depths, extinction must
    correspond to that rescaled coordinate. Scalar DISORT's internal delta-M
    treatment returns physical DFDT and requires physical extinction.
    """
    extinction = _finite(extinction, "extinction")
    if np.any(extinction < 0):
        raise ValueError("extinction must be nonnegative")
    return (
        _finite(dfdt, "dfdt")
        * extinction
        / (_positive(density, "density") * _positive(heat_capacity, "heat_capacity"))
    )


def from_disort(flux, extinction, density, heat_capacity, *, weights=None):
    """Return lower-layer-boundary heating [K/s] from an ARTS ``FluxComponentsProfile``.

    ``extinction`` [1/m] broadcasts to ``flux.dfdt`` (frequency, layer),
    while density and heat capacity normally have shape (layer,).
    Use physical, unscaled extinction for delta-M calculations. Evaluate
    thermodynamic quantities at ``flux.alt_grid[1:]``. Extinction is the
    value within the layer immediately above each output boundary.
    Frequency-dependent extinction is applied *before* spectral integration.
    Optional quadrature ``weights`` are in Hz.
    """
    dfdt = _finite(flux.dfdt, "flux.dfdt")
    if dfdt.shape != (len(flux.freq_grid), len(flux.alt_grid) - 1):
        raise ValueError("flux.dfdt must have shape (frequency, altitude levels - 1)")
    heating = from_optical_depth_derivative(dfdt, extinction, density, heat_capacity)
    if heating.shape != dfdt.shape:
        raise ValueError("thermodynamic inputs must broadcast to flux.dfdt's shape")
    return integrate_spectral(heating, flux.freq_grid, weights=weights)
