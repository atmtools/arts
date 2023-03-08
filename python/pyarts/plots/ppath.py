import pyarts
import numpy as np
import matplotlib.pyplot as plt


def polar_ppath_helper(rad, tht, planetary_radius, rscale, ax=None):
    """ Just the polar plots required by polar_ppath.  This is a helper func

    Parameters
    ----------
    rad : np.array
        List of radiuses [in meters].
    tht : np.array
        List of angles [in radian].
    planetary_radius : float
        A planetary radius in same unit as rad.
    rscale : A rescaler for the radius
        This will rescale values.  If rad is in meters, rscale=1000 means that
        the scale is now in kilometers
    ax : matplotlib axis, optional
        The axis to draw at. The default is None, which generates a default
        polar coordinate.

    Returns
    -------
    ax : As input
        As input.

    """
    if ax is None:
        ax = plt.subplot(111, polar=True)

    ax.plot(tht, rad / rscale + planetary_radius/rscale)
    ax.set_rmin(planetary_radius/rscale)

    ax.set_theta_zero_location('E')
    ax.set_thetalim(-np.pi, np.pi)
    ax.set_thetagrids(np.arange(-180, 179, 30))

    return ax


def polar_ppath_rad_unit(rscale):
    """ Returns the radial unit

    Parameters
    ----------
    rscale : float
        1.0 for "m", 1000 for "km", 1e6 for "Mm.
        Defaults to "???"

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


def polar_ppath_lat(rad, lat, planetary_radius, rscale, ax=None):
    """ Basic plot for ppath latitudes

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
    ax : matplotlib axis, optional
        The axis to draw at. The default is None, which generates a
        default polar coordinate.

    Returns
    -------
    ax : As input.
        As input
    """

    if ax is None:
        ax = plt.subplot(111, polar=True)

    ax = polar_ppath_helper(rad, lat, planetary_radius, rscale, ax)
    ax.set_frame_on(False)
    ax.set_title("Latitude vs "
                 f"{'Altitude' if planetary_radius==0 else 'Radius'}")
    ax.set_thetalim(-np.pi/2, np.pi/2)
    ax.set_thetagrids(np.arange(-90, 91, 45))

    return ax


def polar_ppath_lon(rad, lon, planetary_radius, rscale, ax=None):
    """ Basic plot for ppath longitudes

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
    ax : matplotlib axis, optional
        The axis to draw at. The default is None, which generates a
        default polar coordinate.

    Returns
    -------
    ax : As input.
        As input
    """

    if ax is None:
        ax = plt.subplot(111, polar=True)

    ax = polar_ppath_helper(rad, lon, planetary_radius, rscale, ax)
    ax.set_frame_on(False)
    ax.set_title("Longitude vs "
                 f"{'Altitude' if planetary_radius==0 else 'Radius'}")
    ax.set_theta_zero_location('S')
    ax.set_thetagrids(np.arange(-180, 179, 45))
    ax.set_yticklabels([])  # Disable r-ticks by bad name

    return ax


def polar_ppath_map(lat, lon, planetary_radius, rscale, ax=None):
    """ Basic plot for ppath longitudes

    Parameters
    ----------
    lat : np.array
        List of latitudes [in degrees]
    lon : np.array
        List of longitudes [in degrees].
    planetary_radius : float
        A planetary radius in same unit as rad.
    rscale : A rescaler for the radius
        This will rescale values.  If rad is in meters, rscale=1000
        means that the scale is now in kilometers
    ax : matplotlib axis, optional
        The axis to draw at. The default is None, which generates a
        default polar coordinate.

    Returns
    -------
    ax : As input.
        As input
    """

    if ax is None:
        ax = plt.subplot(111, polar=False)

    plot_data = None
    for [londeg, latdeg] in unwrap_lon(lon, lat):
        if plot_data is None:
            plot_data, = ax.plot(londeg, latdeg)
        else:
            plot_data, = ax.plot(londeg, latdeg, color=plot_data.get_color())

    ax.set_ylim(-90, 90)
    ax.set_xlim(-180, 180)
    ax.set_title("Latitude vs Longitude")
    ax.set_ylabel("Latitude [deg]")
    ax.set_xlabel("Longitude [deg]")
    ax.set_xticks(np.linspace(-180, 180, 7))
    ax.set_yticks(np.linspace(-90, 90, 7))

    return ax


def polar_ppath_za(za, ax=None):
    """ Basic plot for ppath longitudes

    Parameters
    ----------
    za : np.array
        List of Zenith angles [in radians]
    ax : matplotlib axis, optional
        The axis to draw at. The default is None, which generates a
        default polar coordinate.

    Returns
    -------
    ax : As input.
        As input
    """

    if ax is None:
        ax = plt.subplot(111, polar=True)

    ax = polar_ppath_helper(np.ones_like(za), za, 0.0, 1.0, ax)
    ax.set_frame_on(False)
    ax.set_title("Zenith Angle")
    ax.set_thetalim(0, np.pi)
    ax.set_thetagrids(np.arange(0, 181, 45))
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)

    return ax


def polar_ppath_aa(aa, ax=None):
    """ Basic plot for ppath longitudes

    Parameters
    ----------
    aa : np.array
        List of Azimuth angles [in radians]
    ax : matplotlib axis, optional
        The axis to draw at. The default is None, which generates a
        default polar coordinate.

    Returns
    -------
    ax : As input.
        As input
    """

    if ax is None:
        ax = plt.subplot(111, polar=True)

    ax = polar_ppath_helper(np.ones_like(aa), aa, 0.0, 1.0, ax)
    ax.set_frame_on(False)
    ax.set_title("Azimuth Angle")
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.arange(-180, 179, 45))
    ax.set_yticklabels([])  # Disable r-ticks by bad name

    return ax


def unwrap_lon(lon, lat):
    """ Unwraps the lat and lon for plotting purposes.  This is a helper func

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


def polar_ppath(ppath, planetary_radius=0.0, rscale=1000,
                figure_kwargs={"dpi": 300},
                draw_lat_lon=True, draw_map=True, draw_za_aa=False,
                fig=None, axes=None):
    """  Plots a single observation in a polar coordinate system

    The default layout is:
        LAT.vs.RAD | LON.vs.RAD
        -----------------------
              LAT.vs.LON

    The polar plots' central point is at the surface of the planet, i.e., at
    planetary_radius/rscale.  The radius of these plots are the scaled down
    radiuses of the input ppath.pos[0, :] / rscale + planetary_radius/rscale.
    The default radius value is thus just the altitude in kilometers.  If you
    put, e.g., 6371e3 as the planetary_radius, the radius values will be the
    radius from the surface to the highest altitude

    Note also that longitudes are unwrapped, e.g. a step longer than 180
    degrees between Ppath points will wrap around, or rather, create separate
    entries of the lat-lons.

    Parameters
    ----------
    ppath : pyarts.arts.Ppath
        A single propagation path object
    planetary_radius : float, optional
        See polar_part
    rscale : TYPE, optional
        See polar_part
    figure_kwargs : dict, optional
        Arguments to put into plt.figure(). The default is {"dpi": 300}.
    draw_lat_lon : bool, optional
        Whether or not latitude and longitude vs radius angles are drawn.
        Def: True
    draw_map : bool, optional
        Whether or not latitude and longitude map is drawn.  Def: True
    draw_za_aa : bool, optional
        Whether or not Zenith and Azimuth angles are drawn.  Def: False
    fig : A matplotlib figure, optional
        A figure. The default is None, which generates a new figure.
    axes : A list of five matplotlib axes objects, optional
        A tuple of five axis. The default is None, which generates new axes.
        The order are [lat, lon, map, za, aa]

    Returns
    -------
    fig : As input
        As input.
    axes : As input
        As input.

    """
    planetary_radius = float(planetary_radius)

    if fig is None:
        fig = plt.figure(**figure_kwargs)

    if axes is None:
        R = draw_map + (draw_za_aa or draw_lat_lon)
        Z = 2 * draw_lat_lon
        C = 2 * draw_za_aa + Z

        ax_lat = fig.add_subplot(R, C, 1, polar=True) \
            if draw_lat_lon else None
        ax_lon = fig.add_subplot(R, C, 2, polar=True) \
            if draw_lat_lon else None
        ax_za = fig.add_subplot(R, C, 1 + Z, polar=True) \
            if draw_za_aa else None
        ax_aa = fig.add_subplot(R, C, 2 + Z, polar=True) \
            if draw_za_aa else None
        ax_map = fig.add_subplot(R, 1, R, polar=False, aspect=0.5) \
            if draw_map else None

        if draw_za_aa and draw_lat_lon:
            ax_lat.set_position([0.0, 1.0, 0.2, 0.2])
            ax_lon.set_position([0.3, 1.0, 0.2, 0.2])
            ax_za .set_position([0.6, 1.0, 0.2, 0.2])
            ax_aa .set_position([0.9, 1.0, 0.2, 0.2])
            if draw_map:
                ax_map.set_position([0.1, 0.4, 1.0, 0.5])

        axes = [ax_lat, ax_lon, ax_map, ax_za, ax_aa]

    # Set radius and convert degrees
    rad = ppath.pos[:, 0] if ppath.pos.shape[1] > 0 else []
    latdeg = ppath.pos[:, 1] if ppath.pos.shape[1] > 1 else []
    londeg = ppath.pos[:, 2] if ppath.pos.shape[1] > 2 else []
    lat = np.deg2rad(latdeg)
    lon = np.deg2rad(londeg)
    za = np.deg2rad(ppath.los[:, 0]) if ppath.los.shape[1] > 0 else []
    aa = np.deg2rad(ppath.los[:, 1]) if ppath.los.shape[1] > 1 else []

    if draw_lat_lon:
        axes[0] = polar_ppath_lat(rad, lat, planetary_radius, rscale, axes[0])
        axes[0].set_ylabel(f"{'Altitude' if planetary_radius==0 else 'Radius'}"
                           f" [{polar_ppath_rad_unit(rscale)}]")
        axes[1] = polar_ppath_lon(rad, lon, planetary_radius, rscale, axes[1])

    if draw_map:
        axes[2] = polar_ppath_map(latdeg, londeg, planetary_radius, rscale,
                                  axes[2])

    if draw_za_aa:
        axes[3] = polar_ppath_za(za, axes[3])
        axes[3].set_ylabel("Arbitrary unit [-]")
        axes[4] = polar_ppath_aa(aa, axes[4])

    # Return incase people want to modify more
    return fig, axes


def polar_ppath_list(
    ppaths,
    planetary_radius=0.0,
    rscale=1000,
    figure_kwargs={"dpi": 300},
    fig=None,
    axes=None,
    select="end",
    show="poslos"
):
    """ Wraps polar_ppath for a list of ppaths with optional outputs

    This function takes several ppath objects in a list and manipulates them
    based on the option input to form a new ppath object that only has a valid
    pos field.  This new ppath object is the passed directly to polar_ppath
    to plot the polar coordinate and mapping information about pathing

    For example, by default, the option argument is end_pos, which means that
    this call would plot the satellite orbit if the ppaths describe a set of
    satallite (ordered) observations

    Parameters
    ----------
    ppaths : list of Ppath
        A list of path calculations.
    planetary_radius : float, optional
        See polar_part
    rscale : float, optional
        See polar_part
    figure_kwargs : TYPE, optional
        See polar_ppath
    fig : TYPE, optional
        See polar_ppath
    axes : TYPE, optional
        See polar_ppath
    select : str, optional
        The selection criteria for the positions and line of sights in the
        ppath list.  Default is "end".  Options are: "end" for end_pos and
        end_los, "start" for start_pos and start_los, "low" for the lowest
        r pos and los, and "all" for all pos and los
    show : str/list, optional
        Selects what to show.  Default is "poslos" for showing everything.
        The python "in" operator is used to determine what to show, and the
        triggers are: "pos" to show ppath pos, "los" to show ppath los, and
        "no_map" to not show the map.  Note that for all intents and purposes,
        "poslos", "lospos", ["pos", "los"], and "lost purpose" will give the
        exact same result because they all contain "los" and "pos" in a way
        that "in" will find

    Returns
    -------
    See polar_ppath.

    """
    draw_map = "no_map" not in show
    draw_lat_lon = "pos" in show
    draw_za_aa = "los" in show

    my_path = pyarts.arts.Ppath()
    if "end" == select:
        my_path.pos = [ppath.end_pos for ppath in ppaths]
        my_path.los = [ppath.end_los for ppath in ppaths]
    elif "start" == select:
        my_path.pos = [ppath.start_pos for ppath in ppaths]
        my_path.los = [ppath.start_los for ppath in ppaths]
    elif "low" == select:
        my_path.pos = np.concatenate([ppath.pos[ppath.r[:].min() == ppath.r[:]]
                                      for ppath in ppaths])
        my_path.los = np.concatenate([ppath.los[ppath.r[:].min() == ppath.r[:]]
                                      for ppath in ppaths])
    elif "all" == select:
        my_path.pos = np.concatenate([ppath.pos for ppath in ppaths])
        my_path.los = np.concatenate([ppath.los for ppath in ppaths])
    else:
        assert False, f"Unknown selection: '{select}'"

    return polar_ppath(
        ppath=my_path, planetary_radius=planetary_radius, rscale=rscale,
        figure_kwargs=figure_kwargs,
        draw_lat_lon=draw_lat_lon, draw_map=draw_map, draw_za_aa=draw_za_aa,
        fig=fig, axes=axes
    )
