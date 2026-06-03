"""Interactive comparison of the production ``ecef2geodetic_los`` (new)
against the legacy iterative path (old).

A slider drives a geodetic position (latitude, longitude, altitude) and
a zenith angle.  For each azimuth on a fixed grid spanning ``[-180, 180]``
the script round-trips the LOS through ``geodetic_los2ecef`` and then
back through both algorithms, and plots five rows:

* azimuth (input vs new vs old, and ``new - input`` vs ``old - input``),
* zenith angle (same layout),
* recovered altitude (same layout, with input overlaid as a flat line),
* recovered latitude (same layout, with input overlaid as a flat line),
* recovered longitude (same layout, with input overlaid as a flat line).

The "old" path is a Python re-implementation of the iterative
``ecef2geodetic_old`` solver (now retired from the C++ core) plus the
ENH-to-LOS conversion that ``test_geodetic_aa.cpp`` still uses for
comparison.  The "new" path delegates to
``pyarts3.arts.geodetic.ecef2geodetic_los``.
"""

from math import acos, atan, atan2, cos, fabs, hypot, sin, sqrt

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.widgets import Slider

from pyarts3 import arts

DEG2RAD = np.pi / 180.0
RAD2DEG = 180.0 / np.pi

POLELATZZZ = 90.0 - 1e-8
NEAR_POLE_THRESHOLD_ECEF = 1e-15
MAX_ECEF2GEODETIC_ITER = 100


def ecef2geodetic_old(ecef, refellipsoid):
    """Python port of the retired iterative ECEF-to-geodetic solver
    preserved in ``test_geodetic_aa.cpp`` for diagnostic comparison."""
    a, b = refellipsoid
    pos = [0.0, 0.0, 0.0]

    if a == b:
        pos[0] = hypot(ecef[0], ecef[1], ecef[2]) - a
        pos[1] = RAD2DEG * asin_safe(ecef[2] / hypot(ecef[0], ecef[1], ecef[2]))
        pos[2] = RAD2DEG * atan2(ecef[1], ecef[0])
        return pos

    sq = hypot(ecef[0], ecef[1])

    if sq < NEAR_POLE_THRESHOLD_ECEF * a:
        pos[0] = fabs(ecef[2]) - b
        pos[1] = 90.0 if ecef[2] >= 0 else -90.0
        pos[2] = RAD2DEG * atan2(ecef[1], ecef[0])
        return pos

    pos[2] = RAD2DEG * atan2(ecef[1], ecef[0])

    B0 = atan2(ecef[2], sq)
    B = B0 - 1
    N = 0.0
    e2 = 1.0 - (b * b) / (a * a)
    num_iter = 0
    pos[0] = 0.0
    while fabs(B - B0) > 1e-15 and num_iter < MAX_ECEF2GEODETIC_ITER:
        N = a / sqrt(1.0 - e2 * sin(B0) * sin(B0))
        pos[0] = sq / cos(B0) - N
        B = B0
        B0 = atan((ecef[2] / sq) * 1.0 / (1.0 - e2 * N / (N + pos[0])))
        num_iter += 1

    if num_iter == MAX_ECEF2GEODETIC_ITER:
        raise RuntimeError(
            "ECEF to geodetic conversion did not converge. "
            "Input may be too close to singular values."
        )

    pos[1] = RAD2DEG * B0
    return pos


def asin_safe(x):
    return np.arcsin(max(-1.0, min(1.0, x)))


def ecef2geodetic_los_old(ecef, decef, refellipsoid):
    """Legacy LOS round-trip: position via ``ecef2geodetic_old``,
    LOS via the ENU-based formula with the pole azimuth fix."""
    pos = ecef2geodetic_old(ecef, refellipsoid)
    latrad = DEG2RAD * pos[1]
    lonrad = DEG2RAD * pos[2]
    coslat = cos(latrad)
    sinlat = sin(latrad)
    coslon = cos(lonrad)
    sinlon = sin(lonrad)

    enu = [
        -sinlon * decef[0] + coslon * decef[1],
        -sinlat * coslon * decef[0] - sinlat * sinlon * decef[1] + coslat * decef[2],
        coslat * coslon * decef[0] + coslat * sinlon * decef[1] + sinlat * decef[2],
    ]

    twonorm = sqrt(np.dot(enu, enu))
    za = RAD2DEG * acos(enu[2] / twonorm)
    aa = RAD2DEG * atan2(enu[0], enu[1])
    if fabs(pos[1]) > POLELATZZZ:
        aa = RAD2DEG * atan2(decef[1], decef[0])

    return pos, [za, aa]


def ecef2geodetic_los_new(ecef, decef, refellipsoid):
    """Thin wrapper around the production
    ``pyarts3.arts.geodetic.ecef2geodetic_los`` C++ function."""
    pos, los = arts.geodetic.ecef2geodetic_los(
        arts.Vector3(ecef), arts.Vector3(decef), arts.Vector2(refellipsoid)
    )
    return list(pos), list(los)


def azimuth_difference_deg(aa, aa_ref):
    """Wrap-180 aware angular difference in degrees."""
    return ((aa - aa_ref + 180.0) % 360.0) - 180.0


def main():
    wgs84 = arts.planets.Earth.ellipsoid

    fig, axes = plt.subplots(5, 2, figsize=(14, 16), sharex=True)
    ax, ax_diff = axes[0]
    ax_zen, ax_zen_diff = axes[1]
    ax_alt, ax_alt_diff = axes[2]
    ax_lat, ax_lat_diff = axes[3]
    ax_lon, ax_lon_diff = axes[4]
    plt.subplots_adjust(left=0.08, bottom=0.08, right=0.52, hspace=0.4, wspace=0.25)

    aas = np.linspace(-180.0, 180.0, 361 * 5)

    lat0, lon0, alt0, za0 = 0.0, 0.0, 0.0, 90.0
    aa_in = aas.copy()
    za_in = np.full_like(aas, za0)
    alt_in = np.full_like(aas, max(alt0, 1.0))
    lat_in = np.full_like(aas, lat0)
    lon_in = np.full_like(aas, lon0)

    (line_in,) = ax.plot(aas, aa_in, label="Input azimuth", linewidth=2)
    (line_new,) = ax.plot(
        aas, aa_in, "--", label="ecef2geodetic_los (new)", linewidth=2, color="C1"
    )
    (line_old,) = ax.plot(
        aas, aa_in, ":", label="ecef2geodetic_los_old", linewidth=2, color="C2"
    )

    ax.set_xlim(-180.0, 180.0)
    ax.set_ylim(-200.0, 200.0)
    ax.set_ylabel("Azimuth angle [deg]")
    ax.set_title("Round-trip azimuth: input vs new vs old")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    (line_aa_new_diff,) = ax_diff.plot(
        aas, np.zeros_like(aas), "--", label="new - input", linewidth=2, color="C1"
    )
    (line_aa_old_diff,) = ax_diff.plot(
        aas, np.zeros_like(aas), ":", label="old - input", linewidth=2, color="C2"
    )

    ax_diff.set_xlim(-180.0, 180.0)
    ax_diff.set_ylabel("Azimuth diff [deg]")
    ax_diff.set_title("Round-trip azimuth differences")
    ax_diff.grid(True, alpha=0.3)
    ax_diff.legend(loc="upper right")

    (line_za_in,) = ax_zen.plot(aas, za_in, label="Input zenith", linewidth=2)
    (line_za_new,) = ax_zen.plot(
        aas, za_in, "--", label="ecef2geodetic_los (new)", linewidth=2, color="C1"
    )
    (line_za_old,) = ax_zen.plot(
        aas, za_in, ":", label="ecef2geodetic_los_old", linewidth=2, color="C2"
    )

    ax_zen.set_xlim(-180.0, 180.0)
    ax_zen.set_ylim(0.0, 180.0)
    ax_zen.set_ylabel("Zenith angle [deg]")
    ax_zen.set_title("Round-trip zenith: input vs new vs old")
    ax_zen.grid(True, alpha=0.3)

    (line_za_new_diff,) = ax_zen_diff.plot(
        aas, np.zeros_like(aas), "--", label="new - input", linewidth=2, color="C1"
    )
    (line_za_old_diff,) = ax_zen_diff.plot(
        aas, np.zeros_like(aas), ":", label="old - input", linewidth=2, color="C2"
    )

    ax_zen_diff.set_xlim(-180.0, 180.0)
    ax_zen_diff.set_ylabel("Zenith diff [deg]")
    ax_zen_diff.set_title("Round-trip zenith differences")
    ax_zen_diff.grid(True, alpha=0.3)

    (line_alt_in,) = ax_alt.plot(aas, alt_in, label="Input altitude", linewidth=2)
    (line_alt_new,) = ax_alt.plot(
        aas, alt_in, "--", label="ecef2geodetic_los (new)", linewidth=2, color="C1"
    )
    (line_alt_old,) = ax_alt.plot(
        aas, alt_in, ":", label="ecef2geodetic_los_old", linewidth=2, color="C2"
    )

    ax_alt.set_xlim(-180.0, 180.0)
    ax_alt.set_ylabel("Altitude [m]")
    ax_alt.set_title("Round-trip altitude: input vs new vs old")
    ax_alt.grid(True, alpha=0.3)
    ax_alt.set_yscale("linear")

    (line_alt_new_diff,) = ax_alt_diff.plot(
        aas, np.zeros_like(aas), "--", label="new - input", linewidth=2, color="C1"
    )
    (line_alt_old_diff,) = ax_alt_diff.plot(
        aas, np.zeros_like(aas), ":", label="old - input", linewidth=2, color="C2"
    )

    ax_alt_diff.set_xlim(-180.0, 180.0)
    ax_alt_diff.set_ylabel("Altitude diff [m]")
    ax_alt_diff.set_title("Round-trip altitude differences")
    ax_alt_diff.grid(True, alpha=0.3)
    ax_alt_diff.set_yscale("linear")

    (line_lat_in,) = ax_lat.plot(aas, lat_in, label="Input latitude", linewidth=2)
    (line_lat_new,) = ax_lat.plot(
        aas, lat_in, "--", label="ecef2geodetic_los (new)", linewidth=2, color="C1"
    )
    (line_lat_old,) = ax_lat.plot(
        aas, lat_in, ":", label="ecef2geodetic_los_old", linewidth=2, color="C2"
    )

    ax_lat.set_xlim(-180.0, 180.0)
    ax_lat.set_ylim(-90.0, 90.0)
    ax_lat.set_ylabel("Latitude [deg]")
    ax_lat.set_title("Round-trip latitude: input vs new vs old")
    ax_lat.grid(True, alpha=0.3)

    (line_lat_new_diff,) = ax_lat_diff.plot(
        aas, np.zeros_like(aas), "--", label="new - input", linewidth=2, color="C1"
    )
    (line_lat_old_diff,) = ax_lat_diff.plot(
        aas, np.zeros_like(aas), ":", label="old - input", linewidth=2, color="C2"
    )

    ax_lat_diff.set_xlim(-180.0, 180.0)
    ax_lat_diff.set_ylabel("Latitude diff [deg]")
    ax_lat_diff.set_title("Round-trip latitude differences")
    ax_lat_diff.grid(True, alpha=0.3)

    (line_lon_in,) = ax_lon.plot(aas, lon_in, label="Input longitude", linewidth=2)
    (line_lon_new,) = ax_lon.plot(
        aas, lon_in, "--", label="ecef2geodetic_los (new)", linewidth=2, color="C1"
    )
    (line_lon_old,) = ax_lon.plot(
        aas, lon_in, ":", label="ecef2geodetic_los_old", linewidth=2, color="C2"
    )

    ax_lon.set_xlim(-180.0, 180.0)
    ax_lon.set_ylim(-180.0, 180.0)
    ax_lon.set_xlabel("Input azimuth angle [deg]")
    ax_lon.set_ylabel("Longitude [deg]")
    ax_lon.set_title("Round-trip longitude: input vs new vs old")
    ax_lon.grid(True, alpha=0.3)

    (line_lon_new_diff,) = ax_lon_diff.plot(
        aas, np.zeros_like(aas), "--", label="new - input", linewidth=2, color="C1"
    )
    (line_lon_old_diff,) = ax_lon_diff.plot(
        aas, np.zeros_like(aas), ":", label="old - input", linewidth=2, color="C2"
    )

    ax_lon_diff.set_xlim(-180.0, 180.0)
    ax_lon_diff.set_xlabel("Input azimuth angle [deg]")
    ax_lon_diff.set_ylabel("Longitude diff [deg]")
    ax_lon_diff.set_title("Round-trip longitude differences")
    ax_lon_diff.grid(True, alpha=0.3)

    ax_lat_slider = plt.axes((0.63, 0.88, 0.32, 0.025))
    ax_lon_slider = plt.axes((0.63, 0.70, 0.32, 0.025))
    ax_alt_slider = plt.axes((0.63, 0.48, 0.32, 0.025))
    ax_za_slider = plt.axes((0.63, 0.30, 0.32, 0.025))

    ax_alt_slider.set_xscale("log")

    slider_lat = Slider(ax_lat_slider, "Latitude [deg]", -90.0, 90.0, valinit=lat0)
    slider_lon = Slider(ax_lon_slider, "Longitude [deg]", -180.0, 180.0, valinit=lon0)
    slider_alt = Slider(
        ax_alt_slider, "Altitude [m]", 1.0, 1e12, valinit=max(alt0, 1.0), valfmt="%g"
    )
    slider_alt.ax.xaxis.set_major_formatter(
        mticker.FuncFormatter(
            lambda x, _: (
                f"{x:.0f}"
                if x < 1e3
                else f"{x / 1e3:.0f} km"
                if x < 1e6
                else f"{x / 1e6:.0f} Mm"
                if x < 1e9
                else f"{x / 1e9:.0f} Gm"
            )
        )
    )
    slider_za = Slider(ax_za_slider, "Zenith angle [deg]", 0, 180.0, valinit=za0)

    def fmt_altitude(x):
        if abs(x) < 1e3:
            return f"{x:.2f} m"
        if abs(x) < 1e6:
            return f"{x / 1e3:.2f} km"
        if abs(x) < 1e9:
            return f"{x / 1e6:.2f} Mm"
        return f"{x / 1e9:.2f} Gm"

    def update(_):
        lat = slider_lat.val
        lon = slider_lon.val
        alt = slider_alt.val
        za = slider_za.val

        geo = arts.Vector3([alt, lat, lon])
        aa_out_new = np.zeros_like(aas)
        aa_out_old = np.zeros_like(aas)
        za_out_new = np.zeros_like(aas)
        za_out_old = np.zeros_like(aas)
        alt_out_new = np.zeros_like(aas)
        alt_out_old = np.zeros_like(aas)
        lat_out_new = np.zeros_like(aas)
        lat_out_old = np.zeros_like(aas)
        lon_out_new = np.zeros_like(aas)
        lon_out_old = np.zeros_like(aas)

        for i, aa in enumerate(aas):
            los = arts.Vector2([za, aa])
            ecef_los, decef_los = arts.geodetic.geodetic_los2ecef(geo, los, wgs84)
            pos_new, los_new = ecef2geodetic_los_new(
                list(ecef_los), list(decef_los), wgs84
            )
            pos_old, los_old = ecef2geodetic_los_old(
                list(ecef_los), list(decef_los), wgs84
            )
            za_out_new[i] = los_new[0]
            za_out_old[i] = los_old[0]
            aa_out_new[i] = los_new[1]
            aa_out_old[i] = los_old[1]
            alt_out_new[i] = pos_new[0]
            alt_out_old[i] = pos_old[0]
            lat_out_new[i] = pos_new[1]
            lat_out_old[i] = pos_old[1]
            lon_out_new[i] = pos_new[2]
            lon_out_old[i] = pos_old[2]

        line_in.set_ydata(aa_in)
        line_new.set_ydata(aa_out_new)
        line_old.set_ydata(aa_out_old)
        line_za_in.set_ydata(np.full_like(aas, za))
        line_za_new.set_ydata(za_out_new)
        line_za_old.set_ydata(za_out_old)
        line_alt_in.set_ydata(np.full_like(aas, alt))
        line_alt_new.set_ydata(alt_out_new)
        line_alt_old.set_ydata(alt_out_old)
        line_lat_in.set_ydata(np.full_like(aas, lat))
        line_lat_new.set_ydata(lat_out_new)
        line_lat_old.set_ydata(lat_out_old)
        line_lon_in.set_ydata(np.full_like(aas, lon))
        line_lon_new.set_ydata(lon_out_new)
        line_lon_old.set_ydata(lon_out_old)

        aa_new_diff = azimuth_difference_deg(aa_out_new, aas)
        aa_old_diff = azimuth_difference_deg(aa_out_old, aas)
        za_new_diff = za_out_new - za
        za_old_diff = za_out_old - za
        alt_new_diff = alt_out_new - alt
        alt_old_diff = alt_out_old - alt
        lat_new_diff = lat_out_new - lat
        lat_old_diff = lat_out_old - lat
        lon_new_diff = lon_out_new - lon
        lon_old_diff = lon_out_old - lon

        line_aa_new_diff.set_ydata(aa_new_diff)
        line_aa_old_diff.set_ydata(aa_old_diff)
        line_za_new_diff.set_ydata(za_new_diff)
        line_za_old_diff.set_ydata(za_old_diff)
        line_alt_new_diff.set_ydata(alt_new_diff)
        line_alt_old_diff.set_ydata(alt_old_diff)
        line_lat_new_diff.set_ydata(lat_new_diff)
        line_lat_old_diff.set_ydata(lat_old_diff)
        line_lon_new_diff.set_ydata(lon_new_diff)
        line_lon_old_diff.set_ydata(lon_old_diff)

        max_aa = max(
            np.max(np.abs(aa_new_diff)) if aa_new_diff.size else 0.0,
            np.max(np.abs(aa_old_diff)) if aa_old_diff.size else 0.0,
            1e-12,
        )
        ax_diff.set_ylim(-max_aa * 1.1, max_aa * 1.1)

        max_za = max(
            np.max(np.abs(za_new_diff)) if za_new_diff.size else 0.0,
            np.max(np.abs(za_old_diff)) if za_old_diff.size else 0.0,
            1e-12,
        )
        ax_zen_diff.set_ylim(-max_za * 1.1, max_za * 1.1)
        ax_zen_diff.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

        max_alt = max(
            np.max(np.abs(alt_new_diff)) if alt_new_diff.size else 0.0,
            np.max(np.abs(alt_old_diff)) if alt_old_diff.size else 0.0,
            1e-12,
        )
        ax_alt_diff.set_ylim(-max_alt * 1.1, max_alt * 1.1)
        ax_alt_diff.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

        max_lat = max(
            np.max(np.abs(lat_new_diff)) if lat_new_diff.size else 0.0,
            np.max(np.abs(lat_old_diff)) if lat_old_diff.size else 0.0,
            1e-12,
        )
        ax_lat_diff.set_ylim(-max_lat * 1.1, max_lat * 1.1)
        ax_lat_diff.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

        max_lon = max(
            np.max(np.abs(lon_new_diff)) if lon_new_diff.size else 0.0,
            np.max(np.abs(lon_old_diff)) if lon_old_diff.size else 0.0,
            1e-12,
        )
        ax_lon_diff.set_ylim(-max_lon * 1.1, max_lon * 1.1)
        ax_lon_diff.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

        all_alt = np.concatenate([alt_out_new, alt_out_old, [alt]])
        pos_alt = all_alt[all_alt > 0]
        if pos_alt.size:
            lo = max(pos_alt.min(), 1e-3)
            hi = pos_alt.max()
            ax_alt.set_ylim(lo * 0.9, hi * 1.5)

        slider_alt.valtext.set_text(fmt_altitude(alt))
        fig.canvas.draw_idle()

    slider_lat.on_changed(update)
    slider_lon.on_changed(update)
    slider_alt.on_changed(update)
    slider_za.on_changed(update)

    update(None)
    plt.show()


if __name__ == "__main__":
    main()
