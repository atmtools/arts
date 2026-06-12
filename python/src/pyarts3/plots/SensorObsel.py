""" Plotting routine for SensorObsel """

import numpy
import matplotlib
import matplotlib.tri as mtri
import numpy as np
import pyarts3 as pyarts

from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def _line_plot_axes(fig, ax, nkeys):
    return default_fig_ax(fig, ax, 1, nkeys, fig_kwargs={'figsize': (10 * nkeys, 10)})


def _geometry_plot_axes(fig, ax, polar):
    return default_fig_ax(fig,
                          ax,
                          ax_kwargs={"subplot_kw": {'polar': polar}},
                          fig_kwargs={'figsize': (10, 8) if polar else (10, 6)})


def _azimuth_for_plot(azi: np.ndarray, wrap_azimuth: bool) -> np.ndarray:
    if wrap_azimuth:
        return ((azi + 180.0) % 360.0) - 180.0

    return azi


def _los_to_enu(los: np.ndarray) -> np.ndarray:
    zen = np.deg2rad(los[..., 0])
    azi = np.deg2rad(los[..., 1])
    sin_zen = np.sin(zen)
    return np.stack((sin_zen * np.sin(azi),
                     sin_zen * np.cos(azi),
                     np.cos(zen)),
                    axis=-1)


def _antenna_basis(bore_los: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    zen = np.deg2rad(bore_los[0])
    azi = np.deg2rad(bore_los[1])

    cza = np.cos(zen)
    sza = np.sin(zen)
    caa = np.cos(azi)
    saa = np.sin(azi)

    vertical = np.array([-cza * saa, -cza * caa, sza], dtype=float)
    horizontal = np.array([caa, -saa, 0.0], dtype=float)
    bore = np.array([sza * saa, sza * caa, cza], dtype=float)
    return vertical, horizontal, bore


def _local_los_from_enu(enu: np.ndarray, bore_los: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    vertical, horizontal, bore = _antenna_basis(bore_los)

    local_vertical = enu @ vertical
    local_horizontal = enu @ horizontal
    local_bore = np.clip(enu @ bore, -1.0, 1.0)

    zen = np.rad2deg(np.arccos(local_bore))
    azi = np.rad2deg(np.arctan2(local_horizontal, -local_vertical))
    azi = np.where(zen > 0.0, azi, 0.0)
    return zen, azi, local_horizontal, -local_vertical


def _reference_los(data: pyarts.arts.SensorObsel, ifreq: int | None) -> np.ndarray:
    weights = np.asarray(data.weight_matrix, dtype=float)
    iweights = weights[:, :, 0]

    if ifreq is None:
        ref_weights = np.abs(iweights).sum(axis=1)
    else:
        nfreq = iweights.shape[1]
        if ifreq < 0:
            ifreq += nfreq
        ref_weights = np.abs(iweights[:, ifreq])

    los = np.asarray(data.poslos[:, 3:5], dtype=float)
    enu = _los_to_enu(los)

    if np.any(ref_weights > 0.0):
        ref = ref_weights @ enu
    else:
        ref = np.sum(enu, axis=0)

    ref_norm = np.linalg.norm(ref)
    if ref_norm == 0.0:
        return los[0]

    ref /= ref_norm
    ref_zen = np.rad2deg(np.arccos(np.clip(ref[2], -1.0, 1.0)))
    ref_azi = np.rad2deg(np.arctan2(ref[0], ref[1]))
    return np.array([ref_zen, ref_azi], dtype=float)


def _geometry_coordinates(data: pyarts.arts.SensorObsel,
                          frame: str,
                          ifreq: int | None,
                          wrap_azimuth: bool) -> tuple[np.ndarray, np.ndarray, str, str, np.ndarray | None]:
    poslos = np.asarray(data.poslos, dtype=float)
    zen = poslos[:, 3]
    azi = poslos[:, 4]

    if frame == "global":
        plot_azi = _azimuth_for_plot(azi, wrap_azimuth)
        return zen, plot_azi, "Azimuth angle [deg]", "Zenith angle [deg]", None

    if frame != "local":
        raise ValueError("frame must be 'local' or 'global'")

    bore_los = _reference_los(data, ifreq)
    enu = _los_to_enu(poslos[:, 3:5])
    local_zen, local_azi, xoff, yoff = _local_los_from_enu(enu, bore_los)

    radial = local_zen
    theta = _azimuth_for_plot(local_azi, wrap_azimuth)
    _, _, bore = _antenna_basis(bore_los)
    projected_bore = np.clip(enu @ bore, -1.0, 1.0)
    xdeg = np.rad2deg(np.arctan2(xoff, projected_bore))
    ydeg = np.rad2deg(np.arctan2(yoff, projected_bore))

    return radial, theta, "Cross-track offset [deg]", "Along-track offset [deg]", np.stack((xdeg, ydeg), axis=-1)


def _geometry_values(data: pyarts.arts.SensorObsel,
                     pol: pyarts.arts.Stokvec,
                     ifreq: int | None) -> tuple[np.ndarray, str]:
    if ifreq is None:
        values = np.asarray(data.weight_matrix.reduce(pol,
                                                      along_poslos=False,
                                                      along_freq=True),
                            dtype=float)
        return values, "Integrated weight"

    weights = np.asarray(data.weight_matrix, dtype=float)
    nfreq = weights.shape[1]
    if ifreq < -nfreq or ifreq >= nfreq:
        raise IndexError(f"ifreq {ifreq} out of range for {nfreq} frequencies")

    if ifreq < 0:
        ifreq += nfreq

    pol_vec = np.asarray(pol, dtype=float)
    values = weights[:, ifreq, :] @ pol_vec
    return values, f"Weight at f[{ifreq}] = {float(data.f_grid[ifreq]):g}"


def _geometry_plot_type(point_spread: bool, type: str | None) -> str:
    if type is None:
        return "points"

    plot_type = str(type).lower()
    if plot_type == "scatter":
        plot_type = "points"

    valid_types = ("points", "cloud", "contour")
    if plot_type not in valid_types:
        choices = ", ".join(repr(choice) for choice in valid_types)
        raise ValueError(f"type must be one of {choices}")

    if plot_type != "points" and not point_spread:
        raise ValueError(f"type={type!r} requires point_spread=True")

    return plot_type


def _triangulated_point_spread(local_xy: np.ndarray,
                               values: np.ndarray) -> tuple[mtri.Triangulation, np.ndarray]:
    x = np.asarray(local_xy[:, 0], dtype=float)
    y = np.asarray(local_xy[:, 1], dtype=float)
    values = np.asarray(values, dtype=float)

    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(values)
    x = x[finite]
    y = y[finite]
    values = values[finite]

    if x.size < 3:
        raise ValueError("Interpolated point-spread plots require at least 3 geometry samples")

    coords = np.stack((x, y), axis=-1)
    unique_coords, inverse = np.unique(coords, axis=0, return_inverse=True)
    if unique_coords.shape[0] != coords.shape[0]:
        unique_values = np.zeros(unique_coords.shape[0], dtype=float)
        counts = np.zeros(unique_coords.shape[0], dtype=int)
        np.add.at(unique_values, inverse, values)
        np.add.at(counts, inverse, 1)
        x = unique_coords[:, 0]
        y = unique_coords[:, 1]
        values = unique_values / counts

    if x.size < 3:
        raise ValueError("Interpolated point-spread plots require at least 3 unique geometry samples")

    triangulation = mtri.Triangulation(x, y)
    if triangulation.triangles.size == 0:
        raise ValueError("Interpolated point-spread plots require non-collinear geometry samples")

    return triangulation, values


def _contour_levels(values: np.ndarray, levels: int | np.ndarray) -> np.ndarray | int:
    if np.ndim(levels) != 0:
        return levels

    nlevels = max(int(levels), 2)
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    if np.isclose(vmin, vmax):
        delta = 1e-6 if vmin == 0.0 else abs(vmin) * 1e-6
        return np.array((vmin - delta, vmax + delta), dtype=float)

    return np.linspace(vmin, vmax, nlevels)


def _plot_geometry(data: pyarts.arts.SensorObsel,
                   *,
                   fig: matplotlib.figure.Figure | None,
                   ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None,
                   pol: pyarts.arts.Stokvec,
                   polar: bool,
                   point_spread: bool,
                   type: str | None,
                   ifreq: int | None,
                   colorbar: bool,
                   frame: str,
                   wrap_azimuth: bool,
                   **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    fig, ax = _geometry_plot_axes(fig, ax, polar)
    axis = select_flat_ax(ax, 0)

    radial, theta, xlabel, ylabel, local_xy = _geometry_coordinates(data,
                                                                    frame,
                                                                    ifreq,
                                                                    wrap_azimuth)
    plot_type = _geometry_plot_type(point_spread, type)

    if plot_type == "points":
        default_kwargs = {'marker': 'o', 'linestyle': 'None'}
        for key, value in default_kwargs.items():
            kwargs.setdefault(key, value)

    artist = None
    if point_spread:
        values, label = _geometry_values(data, pol, ifreq)
        kwargs.setdefault('cmap', 'viridis')

        if plot_type == "cloud" or plot_type == "contour":
            if polar:
                raise ValueError("Interpolated point-spread types require polar=False")
            if frame != "local" or local_xy is None:
                raise ValueError("Interpolated point-spread types require frame='local'")

            triangulation, values = _triangulated_point_spread(local_xy, values)
            if plot_type == "cloud":
                kwargs.setdefault('shading', 'gouraud')
                artist = axis.tripcolor(triangulation, values, **kwargs)
            else:
                levels = _contour_levels(values, kwargs.pop('levels', 32))
                artist = axis.tricontourf(triangulation, values, levels=levels, **kwargs)
        elif polar:
            artist = axis.scatter(np.deg2rad(theta), radial, c=values, **kwargs)
        elif frame == "local":
            artist = axis.scatter(local_xy[:, 0], local_xy[:, 1], c=values, **kwargs)
        else:
            artist = axis.scatter(theta, radial, c=values, **kwargs)
        if colorbar:
            fig.colorbar(artist, ax=axis, label=label)
    else:
        if polar:
            artist = axis.plot(np.deg2rad(theta), radial, **kwargs)
        elif frame == "local":
            artist = axis.plot(local_xy[:, 0], local_xy[:, 1], **kwargs)
        else:
            artist = axis.plot(theta, radial, **kwargs)

    if polar:
        axis.set_theta_zero_location("N")
        axis.set_theta_direction(-1)
        axis.set_rlabel_position(225)
        axis.set_ylim(0.0, np.max(radial) if radial.size else 1.0)
    else:
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        if frame == "local":
            axis.set_aspect('equal', adjustable='box')

    return fig, ax


def plot(data: pyarts.arts.SensorObsel,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         keys: str | list = "f",
         pol: str | pyarts.arts.Stokvec = "I",
         polar: bool = False,
         point_spread: bool = False,
         type: str | None = None,
         ifreq: int | None = None,
         colorbar: bool = True,
         frame: str = "global",
         wrap_azimuth: bool = False,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot a sensor observational element.

    By default, this follows the same line-plot convention as
    :mod:`pyarts3.plots.ArrayOfSensorObsel`: the selected quantity is reduced
    over the non-plotted axis and drawn against either frequency or one sensor
    geometry component.

    Setting ``polar=True`` switches to a geometry view.  By default, geometry
    plots use the stored ``SensorObsel.poslos`` zenith and azimuth values so
    the plotted coordinates match the observation element directly.
    Setting ``point_spread=True`` colors the geometry samples by their reduced
    weights, integrated over all frequencies unless ``ifreq`` selects one
    frequency column.  The geometry ``type`` selects whether those weighted
    samples appear as discrete points, a smooth triangulated cloud, or filled
    triangulated contours.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        ant = pyarts.arts.sensor.GaussianAiryAntenna(np.linspace(0, 1, 11),
                                                     np.linspace(0, 330, 12),
                                                     0.3)
        ch = pyarts.arts.sensor.DiracChannel(100e9)
        obsel = ant(ch, [1, 0, 0], [90, 0])

        pyarts.plots.SensorObsel.plot(obsel, polar=True, point_spread=True)
        pyarts.plots.SensorObsel.plot(obsel, point_spread=True, wrap_azimuth=True)
        pyarts.plots.SensorObsel.plot(obsel, point_spread=True, frame="local")
        pyarts.plots.SensorObsel.plot(obsel, point_spread=True, frame="local", type="cloud")
        pyarts.plots.SensorObsel.plot(obsel, point_spread=True, frame="local", type="contour")

    Parameters
    ----------
    data : ~pyarts3.arts.SensorObsel
        A single sensor observation element.
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    keys : str | list
        The keys to use for line plotting. Options are in :class:`~pyarts3.arts.SensorKeyType`.
        Ignored when ``polar`` or ``point_spread`` is enabled.
    pol : str | ~pyarts3.arts.Stokvec
        The polarization to use for plotting. Defaults to ``"I"``.
    polar : bool, optional
        If True, plot the observation geometry in polar coordinates. Defaults to False.
    point_spread : bool, optional
        If True, color the observation geometry by weight values. Defaults to False.
    type : str | None, optional
        Geometry rendering type. ``"points"`` shows discrete samples,
        ``"cloud"`` draws a smooth triangulated color cloud, and
        ``"contour"`` draws filled triangulated contours. Interpolated types
        require ``point_spread=True``, ``frame="local"``, and ``polar=False``.
        Defaults to ``None``, which behaves like ``"points"``.
    ifreq : int | None, optional
        Frequency index to use for point-spread values. If None, integrates over frequency.
    colorbar : bool, optional
        If True, add a colorbar to point-spread plots. Defaults to True.
    frame : str, optional
        Geometry frame for ``polar`` and ``point_spread`` plots. ``"global"``
        uses the stored zenith/azimuth coordinates from ``SensorObsel.poslos``,
        while ``"local"`` uses a bore-centered frame inferred from the weighted
        LOS samples. Defaults to ``"global"``.
    wrap_azimuth : bool, optional
        If True, wrap azimuths to ``[-180, 180)`` before geometry plotting.
        Defaults to False.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input. Otherwise the created Figure.
    ax :
        As input if input. Otherwise the created Axes.
    """

    pol = pyarts.arts.Stokvec(pol)

    if type is not None and not (polar or point_spread):
        _geometry_plot_type(point_spread, type)

    if polar or point_spread:
        return _plot_geometry(data,
                              fig=fig,
                              ax=ax,
                              pol=pol,
                              polar=polar,
                              point_spread=point_spread,
                              type=type,
                              ifreq=ifreq,
                              colorbar=colorbar,
                              frame=frame,
                              wrap_azimuth=wrap_azimuth,
                              **kwargs)

    keys = [keys] if isinstance(keys, str) else keys
    nkeys = len(keys)
    if nkeys == 0:
        return fig, ax

    fig, ax = _line_plot_axes(fig, ax, nkeys)

    key_map = {
        pyarts.arts.SensorKeyType.f: None,
        pyarts.arts.SensorKeyType.alt: 0,
        pyarts.arts.SensorKeyType.lat: 1,
        pyarts.arts.SensorKeyType.lon: 2,
        pyarts.arts.SensorKeyType.zen: 3,
        pyarts.arts.SensorKeyType.azi: 4,
    }

    for isub, raw_key in enumerate(keys):
        key = pyarts.arts.SensorKeyType(raw_key)
        idx = key_map[key]

        values = data.weight_matrix.reduce(pol,
                                           along_poslos=idx is None,
                                           along_freq=idx is not None)
        x = data.f_grid if idx is None else data.poslos[:, idx]
        select_flat_ax(ax, isub).plot(x, values, **kwargs)

    return fig, ax