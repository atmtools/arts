import pyarts
import numpy as np
import matplotlib.pyplot as plt


def polar_part(rad, tht, planetary_radius, rscale, ax=None):
    """ Just the polar plots required by polar_ppath.  This is a helper func

    Parameters
    ----------
    rad : np.array
        List of radiuses [in rscale=1 units].
    tht : np.array
        List of angles [in radian].
    planetary_radius : float
        A planetary radius in same unit as rad.
    rscale : A rescaler for the radius
        This will rescale values.  If rad is in meters, rscale=1000 means that
        the scale is now in kilometers
    ax : matplotlib axis, optional
        The axis to draw at. The default is None, which generates a default
        polar coordinat.

    Returns
    -------
    ax : As input
        As input.

    """
    if ax is None:
        ax = plt.subplot(111, projection='polar')

    ax.plot(tht, rad / rscale + planetary_radius/rscale, 'b')
    ax.set_rmin(planetary_radius/rscale)

    ax.set_theta_zero_location('E')
    ax.set_thetalim(-np.pi, np.pi)
    ax.set_thetagrids(np.arange(-180, 179, 45))

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
                figure_kwargs={"dpi": 300}, fig=None, axes=None):
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
    fig : A matplotlib figure, optional
        A figure. The default is None, which generates a new figure.
    axes : A list of three matplotlib axes objects, optional
        A tuple of three axis. The default is None, which generates new axes.

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
        axes = [fig.add_subplot(2, 2, 1, polar=True),
                fig.add_subplot(2, 2, 2, polar=True),
                fig.add_subplot(2, 1, 2, polar=False)]

    # Set radius and convert degrees
    rad = ppath.pos[:, 0]
    latdeg = ppath.pos[:, 1]
    londeg = ppath.pos[:, 2]
    lat = np.deg2rad(ppath.pos[:, 1])
    lon = np.deg2rad(ppath.pos[:, 2])

    # Plot data for LAT.v.R
    ax1 = polar_part(rad, lat, planetary_radius, rscale, axes[0])
    ax1.set_frame_on(False)
    ax1.set_title("Latitude vs Radius")
    ax1.set_thetalim(-np.pi/2, np.pi/2)
    ax1.set_thetagrids(np.arange(-90, 91, 45))

    # Plot data for LON.v.R
    ax2 = polar_part(rad, lon, planetary_radius, rscale, axes[1])
    ax2.set_frame_on(False)
    ax2.set_title("Longitude vs Radius")
    ax2.set_yticklabels([])  # Disable r-ticks by bad name

    # Plot data for LAT.v.LON
    ax3 = axes[2]
    for [londeg, latdeg] in unwrap_lon(londeg, latdeg):
        ax3.plot(londeg, latdeg, 'b')
    ax3.set_ylim(-90, 90)
    ax3.set_xlim(-180, 180)
    ax3.set_title("Latitude vs Longitude")
    ax3.set_ylabel("Latitude [deg]")
    ax3.set_xlabel("Longitude [deg]")

    # Return incase people want to modify more
    return fig, [ax1, ax2, ax3]


def polar_ppath_list(
    ppaths,
    planetary_radius=0.0,
    rscale=1000,
    figure_kwargs={"dpi": 300},
    fig=None,
    axes=None,
    option="end_pos",
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
    option : str, optional
        The options changes how the plotting happens. The default is "end_pos".

        Accepted input values and the actions are:

            end_pos : use end_pos from all ppath objects

            start_pos : use start_pos from all ppath objects

            pos : Concatenate all pos from all ppath objects

            lowest : Concatenate all lowest radius pos from all ppath objects

    Returns
    -------
    See polar_ppath.

    """
    my_path = pyarts.arts.Ppath()
    if option.lower() == "end_pos":
        my_path.pos = [ppath.end_pos for ppath in ppaths]
    elif option.lower() == "start_pos":
        my_path.pos = [ppath.start_pos for ppath in ppaths]
    elif option.lower() == "pos":
        my_path.pos = np.concatenate([ppath.pos for ppath in ppaths])
    elif option.lower() == "lowest":
        my_path.pos = np.concatenate([ppath.pos[ppath.r[:].min() == ppath.r[:]]
                                      for ppath in ppaths])
    else:
        assert False, f"Not a valid option: '{option}'"

    return polar_ppath(
        my_path, planetary_radius, rscale, figure_kwargs, fig, axes
    )
