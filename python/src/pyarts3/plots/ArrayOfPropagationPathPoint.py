""" Plotting routine the propagation path in polar coordinates """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def polar_ray_path_helper(rad, tht, planetary_radius, rscale, ax=None):
    """Just the polar plots required by :func:`polar_ray_path`.

    Parameters
    ----------
    rad : np.array
        List of radiuses [in meters].
    tht : np.array
        List of angles [in radian].
    planetary_radius : float
        A planetary radius in same unit as rad.
    rscale : float
        This will rescale values.  As rad is in meters, rscale=1000 means that
        the scale is now in kilometers.  See :func:`polar_ray_path_rad_unit` for
        good options
    ax : subs, optional
        The axis to draw at. The default is None, which generates a default
        polar coordinate.

    Returns
    -------
    ax : As input
        As input.

    """
    if ax is None:
        ax = plt.subplot(111, polar=True)

    st = '-' if len(rad) > 1 else 'x'

    ax.plot(tht, rad / rscale + planetary_radius / rscale, st)
    ax.set_rmin(planetary_radius / rscale)

    ax.set_theta_zero_location("E")
    ax.set_thetalim(-np.pi, np.pi)
    ax.set_thetagrids(np.arange(-180, 179, 30))

    return ax


def polar_ray_path_rad_unit(rscale):
    """Returns the radial unit

    Parameters
    ----------
    rscale : float
        1.0 for "m", 1000 for "km", 1e6 for "Mm.
        Otherwise returns "???"

    Returns
    -------
    str
        The unit or "???".

    """
    if rscale == 1.0:
        return "m"
    elif rscale == 1000.0:
        return "km"
    elif rscale == 1e6:
        return "Mm"
    else:
        return "???"


def polar_ray_path_lat(rad, lat, planetary_radius, rscale, ax=None):
    """Basic plot for ray_path latitudes

    Parameters
    ----------
    rad : np.array
        List of radiuses [in meters]
    lat : np.array
        List of latitudes [in radian].
    planetary_radius : float
        A planetary radius in same unit as rad.
    rscale : A rescaler for the radius
        This will rescale values.  If rad is in meters, rscale=1000
        means that the scale is now in kilometers
    ax : subs, optional
        The axis to draw at. The default is None, which generates a
        default polar coordinate.

    Returns
    -------
    ax : As input.
        As input
    """

    if ax is None:
        ax = plt.subplot(111, polar=True)

    ax = polar_ray_path_helper(rad, lat, planetary_radius, rscale, ax)
    ax.set_frame_on(False)
    ax.set_title(
        "Latitude vs " f"{'Altitude' if planetary_radius == 0 else 'Radius'}"
    )
    ax.set_thetalim(-np.pi / 2, np.pi / 2)
    ax.set_thetagrids(np.arange(-90, 91, 45))

    return ax


def polar_ray_path_lon(rad, lon, planetary_radius, rscale, ax=None):
    """Basic plot for ray_path longitudes

    Parameters
    ----------
    rad : np.array
        List of radiuses [in meters]
    lon : np.array
        List of longitudes [in radian].
    planetary_radius : float
        A planetary radius in same unit as rad.
    rscale : A rescaler for the radius
        This will rescale values.  If rad is in meters, rscale=1000
        means that the scale is now in kilometers
    ax : subs, optional
        The axis to draw at. The default is None, which generates a
        default polar coordinate.

    Returns
    -------
    ax : As input.
        As input
    """

    if ax is None:
        ax = plt.subplot(111, polar=True)

    ax = polar_ray_path_helper(rad, lon, planetary_radius, rscale, ax)
    ax.set_frame_on(False)
    ax.set_title(
        "Longitude vs " f"{'Altitude' if planetary_radius == 0 else 'Radius'}"
    )
    ax.set_theta_zero_location("S")
    ax.set_thetagrids(np.arange(-180, 179, 45))
    ax.set_yticklabels([])  # Disable r-ticks by bad name

    return ax


def polar_ray_path_map(lat, lon, ax=None):
    """Basic plot for ray_path longitudes

    Parameters
    ----------
    lat : np.array
        List of latitudes [in degrees]
    lon : np.array
        List of longitudes [in degrees].
    ax : subs, optional
        The axis to draw at. The default is None, which generates a
        default polar coordinate.

    Returns
    -------
    ax : As input.
        As input
    """

    if ax is None:
        ax = plt.subplot(111, polar=False)

    st = '-' if len(lat) > 1 else 'x'

    plot_data = None
    for [londeg, latdeg] in unwrap_lon(lon, lat):
        if plot_data is None:
            (plot_data,) = ax.plot(londeg, latdeg, st)
        else:
            (plot_data,) = ax.plot(londeg, latdeg, st,
                                   color=plot_data.get_color())

    ax.set_ylim(-90, 90)
    ax.set_xlim(-180, 180)
    ax.set_title("Latitude vs Longitude")
    ax.set_ylabel("Latitude [deg]")
    ax.set_xlabel("Longitude [deg]")
    ax.set_xticks(np.linspace(-180, 180, 7))
    ax.set_yticks(np.linspace(-90, 90, 7))

    return ax


def polar_ray_path_za(za, ax=None):
    """Basic plot for ray_path zeniths

    Parameters
    ----------
    za : np.array
        List of Zenith angles [in radians]
    ax : subs, optional
        The axis to draw at. The default is None, which generates a
        default polar coordinate.

    Returns
    -------
    ax : As input.
        As input
    """

    if ax is None:
        ax = plt.subplot(111, polar=True)

    ax = polar_ray_path_helper(np.ones_like(za), za, 0.0, 1.0, ax)
    ax.set_frame_on(False)
    ax.set_title("Zenith Angle")
    ax.set_thetalim(0, np.pi)
    ax.set_thetagrids(np.arange(0, 181, 45))
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    return ax


def polar_ray_path_aa(aa, ax=None):
    """Basic plot for ray_path azimuths

    Parameters
    ----------
    aa : np.array
        List of Azimuth angles [in radians]
    ax : subs, optional
        The axis to draw at. The default is None, which generates a
        default polar coordinate.

    Returns
    -------
    ax : As input.
        As input
    """

    if ax is None:
        ax = plt.subplot(111, polar=True)

    ax = polar_ray_path_helper(np.ones_like(aa), aa, 0.0, 1.0, ax)
    ax.set_frame_on(False)
    ax.set_title("Azimuth Angle")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.arange(-180, 179, 45))
    ax.set_yticklabels([])  # Disable r-ticks by bad name

    return ax


def unwrap_lon(lon, lat):
    """Unwraps the lat and lon for plotting purposes.  This is a helper func

    Parameters
    ----------
    lon : np.array
        A list of latitudes.
    lat : np.array
        A list of longitudes.

    Returns
    -------
    list
        one or more pairs of [lon, lat] that when plotted does not allow
        wrapping.

    """
    jumps = np.nonzero(np.abs(np.diff(lon)) > 180)[0]

    if len(jumps) == 0:
        return [[lon, lat]]

    out = []
    i = 0
    for ind in jumps:
        ind = ind + 1
        out.append([lon[i:ind], lat[i:ind]])
        i = ind
    out.append([lon[i:-1], lat[i:-1]])
    return out


def polar_ray_path_default_figure(figure_kwargs):
    """Get the default figure by standard inputs

    Parameters
    ----------
    figure_kwargs : dict, optional
        Arguments to put into plt.figure().

    Returns
    -------
    A matplotlib figure, optional
        A figure
    """
    return plt.figure(**figure_kwargs)


def polar_ray_path_default_subs(fig, draw_lat_lon, draw_map, draw_za_aa):
    """Get the default subs

    Parameters
    ----------
    fig : Figure
        A figure
    draw_lat_lon : bool, optional
        Whether or not latitude and longitude vs radius angles are drawn.
        Def: True
    draw_map : bool, optional
        Whether or not latitude and longitude map is drawn.  Def: True
    draw_za_aa : bool, optional
        Whether or not Zenith and Azimuth angles are drawn.  Def: False

    Returns
    -------
    A list of five subs
        A tuple of five axis.

        The order is [lat, lon, map, za, aa]

    """
    R = draw_map + (draw_za_aa or draw_lat_lon)
    Z = 2 * draw_lat_lon
    C = 2 * draw_za_aa + Z

    ax_lat = fig.add_subplot(R, C, 1, polar=True) if draw_lat_lon else None
    ax_lon = fig.add_subplot(R, C, 2, polar=True) if draw_lat_lon else None
    ax_za = fig.add_subplot(R, C, 1 + Z, polar=True) if draw_za_aa else None
    ax_aa = fig.add_subplot(R, C, 2 + Z, polar=True) if draw_za_aa else None
    ax_map = (
        fig.add_subplot(R, 1, R, polar=False, aspect=0.5) if draw_map else None
    )

    if draw_za_aa and draw_lat_lon:
        ax_lat.set_position([0.0, 1.0, 0.2, 0.2])
        ax_lon.set_position([0.3, 1.0, 0.2, 0.2])
        ax_za.set_position([0.6, 1.0, 0.2, 0.2])
        ax_aa.set_position([0.9, 1.0, 0.2, 0.2])
        if draw_map:
            ax_map.set_position([0.1, 0.4, 1.0, 0.5])

    return [ax_lat, ax_lon, ax_map, ax_za, ax_aa]


def plot(
    ray_path: pyarts.arts.ArrayOfPropagationPathPoint,
    *,
    planetary_radius: float = 0.0,
    rscale: float = 1000,
    figure_kwargs: dict = {"dpi": 300},
    draw_lat_lon: bool = True,
    draw_map: bool = True,
    draw_za_aa: bool = False,
    fig=None,
    subs=None,
    **kwargs
):
    """Plots a single observation in a polar coordinate system

    Use the draw_* variables to select which plots are done

    The polar plots' central point is at the surface of the planet, i.e., at
    planetary_radius/rscale.  The radius of these plots are the scaled down
    radiuses of the input ray_path[0].pos[0] / rscale + planetary_radius/rscale.
    The default radius value is thus just the altitude in kilometers.  If you
    put, e.g., 6371e3 as the planetary_radius, the radius values will be the
    radius from the surface to the highest altitude

    Note also that longitudes are unwrapped, e.g. a step longer than 180
    degrees between ray_path points will wrap around, or rather, create separate
    entries of the lat-lons.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        ws = pyarts.Workspace()

        ws.atmospheric_fieldRead(toa=100e3, basename="planets/Earth/afgl/tropical/")
        ws.surface_fieldEarth()
        ws.ray_path_observer_agendaSetGeometric(
            add_crossings=True, remove_non_crossings=True
        )
        ws.ray_path_observersFieldProfilePseudo2D(nup=3, nlimb=3, ndown=3)
        ws.ray_path_fieldFromObserverAgenda()

        f, a = None, None
        for x in ws.ray_path_field:
            f, a = pyarts.plots.ArrayOfPropagationPathPoint.plot(
                x, draw_za_aa=True, draw_map=False, fig=f, subs=a
            )

    Parameters
    ----------
    ray_path : ~pyarts3.arts.ArrayOfPropagationPathPoint
        A single propagation path object
    planetary_radius : float, optional
        See ``polar_ray_path_helper`` in source tree
    rscale : float, optional
        See ``polar_ray_path_helper`` in source tree
    figure_kwargs : dict, optional
        Arguments to put into plt.figure(). The default is {"dpi": 300}.
    draw_lat_lon : bool, optional
        Whether or not latitude and longitude vs radius angles are drawn.
        Def: True
    draw_map : bool, optional
        Whether or not latitude and longitude map is drawn.  Def: True
    draw_za_aa : bool, optional
        Whether or not Zenith and Azimuth angles are drawn.  Def: False
    fig : Figure, optional
        A figure. The default is None, which generates a new figure.
    subs : A list of five subaxis, optional
        A tuple of five subaxis. The default is None, which generates new subs.
        The order is [lat, lon, map, za, aa]

    Returns
    -------
    fig : As input
        As input.
    subs : As input
        As input.

    """
    if fig is None:
        fig = polar_ray_path_default_figure(figure_kwargs)

    if subs is None:
        subs = polar_ray_path_default_subs(
            fig, draw_lat_lon, draw_map, draw_za_aa
        )

    # Set radius and convert degrees
    rad = np.array([x.pos[0] for x in ray_path])
    latdeg = np.array([x.pos[1] for x in ray_path])
    londeg = np.array([x.pos[2] for x in ray_path])
    zadeg = np.array([x.los[0] for x in ray_path])
    aadeg = np.array([x.los[1] for x in ray_path])

    lat = np.deg2rad(latdeg)
    lon = np.deg2rad(londeg)
    za = np.deg2rad(zadeg)
    aa = np.deg2rad(aadeg)

    if draw_lat_lon:
        subs[0] = polar_ray_path_lat(rad, lat, planetary_radius, rscale, subs[0])
        subs[0].set_ylabel(
            f"{'Altitude' if planetary_radius == 0 else 'Radius'}"
            f" [{polar_ray_path_rad_unit(rscale)}]"
        )
        subs[1] = polar_ray_path_lon(rad, lon, planetary_radius, rscale, subs[1])

    if draw_map:
        subs[2] = polar_ray_path_map(latdeg, londeg, subs[2])

    if draw_za_aa:
        subs[3] = polar_ray_path_za(za, subs[3])
        subs[3].set_ylabel("Arbitrary unit [-]")
        subs[4] = polar_ray_path_aa(aa, subs[4])

    # Return incase people want to modify more
    return fig, subs
