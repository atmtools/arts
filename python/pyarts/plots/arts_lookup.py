"""Visualize an absorption lookup table.

Author: oliver.lemke@uni-hamburg.de
"""
import re
from itertools import zip_longest

import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d
from scipy import constants

__all__ = [
    'plot_arts_lookup',
]

_molar_mass_dry_air = 28.9645e-3  # kg mol^-1
_gas_constant_dry_air = constants.gas_constant / _molar_mass_dry_air  # J K^-1 kg^-1


def _calc_lookup_species_count(lookup):
    """Calculate number of cross sections per species.

    Usually one, except for the nonlinear species.
    """
    nlsspecies = [int(s) for s in lookup.nonlinspecs.data]
    speciescount = np.ones_like(np.array(lookup.specs, dtype=object),
                                dtype=int)
    if len(nlsspecies) != 0:
        speciescount[nlsspecies] = lookup.nls_pert.data.size
    return speciescount


def _get_lookup_species_index(lookup, species, vmrpert):
    """Get index of given species in lookup table."""
    ret = 0
    spindex = lookup.specs.data.index(species)
    nlsspecies = lookup.nonlinspecs.data
    speciescount = _calc_lookup_species_count(lookup)
    if nlsspecies is not None and spindex in nlsspecies:
        if vmrpert >= speciescount[spindex]:
            raise RuntimeError(
                'Nonlinear species VMR perturbation index too large')
        ret = vmrpert

    return ret + (np.sum(speciescount[0:spindex]) if spindex > 0 else 0)


def _add_opacity_legend(ax=None):
    """Add legend to an opacity lookup table plot."""
    if ax is None:
        ax = plt.gca()

    blue_line = Line2D([], [], label='species opacity')
    black_line = Line2D([], [], color='k', linewidth=1., label='total opacity')
    dashed_line = Line2D([], [],
                         color='k',
                         linestyle='--',
                         linewidth=1.,
                         label='opacity=1')

    handles = [blue_line, black_line, dashed_line]
    labels = [h.get_label() for h in handles]

    ax.legend(handles=handles,
              labels=labels,
              fontsize='xx-small',
              loc='upper left',
              ncol=6)


def _add_xsec_legend(lookup, ipressures, ax=None):
    """Add legend to a cross section lookup table plot."""
    if ax is None:
        ax = plt.gca()

    pgrid = lookup.p_grid.data
    colors = [plt.cm.viridis(i) for i in np.linspace(0, 1, len(ipressures))]
    handles = [
        Line2D([], [], color=colors[i], label=f'{pgrid[ip]/100.:8.3f} hPa')
        for i, ip in enumerate(ipressures)
    ]

    labels = [h.get_label() for h in handles]

    ax.legend(handles=handles,
              labels=labels,
              fontsize='xx-small',
              loc='upper left',
              ncol=6)


def _setup_lookup_figure(lookup, cols=3, species=None):
    """Create the figure and axes objects for the lookup table plot."""
    if species is None:
        species = lookup.specs
    rows = int(np.ceil(len(species) / cols))
    fig, ax = plt.subplots(rows + 1, cols, figsize=(4 * cols, (rows + 1) * 2))
    fig.tight_layout()

    return rows, cols, fig, ax


def plot_lookup_xsec(lookup,
                     ipressures,
                     species=None,
                     ax=None,
                     tpert=0,
                     vmrpert=0):
    """Plot the cross section for one or more species of an ARTS lookup table.

    Parameters:
        lookup (pyarts.classes.GasAbsLookup): ARTS lookup table.
        ipressures (ndarray(int)): Indices of pressure levels to plot.
        species (pyarts.classes.ArrayOfArrayOfSpeciesTag, optional):
            ARTS species tags. If none is given, plots all species in the lookup
            table for the given vmr perturbation.
        ax (AxesSubplot, optional): Axes to plot in.
        vmrpert (int): Index of vmr perturbation for nonlinear species to plot.
        tpert (int): Index of temperature perturbation to plot.
    """
    if ax is None:
        ax = plt.gca()

    ax.set_yscale('log')
    if species is None:
        species = lookup.specs.data

    for tag in species:
        ax.set_prop_cycle(
            cycler('color', [
                plt.cm.viridis(i) for i in np.linspace(0, 1, len(ipressures))
            ]))
        for pi in ipressures:
            xsec = lookup.xsec.data[
                tpert,
                _get_lookup_species_index(lookup, tag, vmrpert), :, pi]
            ax.plot(lookup.f_grid, xsec, label=f'{pi/100.:8.3f} hPa')

    if len(species) > 1:
        ax.legend(fontsize='xx-small', frameon=False)
    else:
        ax.set_title(',\n'.join(
            re.sub(r'(-\*)+$', '', str(s.as_string)) for s in species[0]),
                     y=1. - len(species[0]) * 0.05,
                     fontsize='xx-small')

    def formatter(x, pos):
        return rf'${x/1e9:g}$'

    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_minor_formatter(formatter)
    ax.tick_params(axis='both', which='major', labelsize='xx-small')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    return ax


def plot_lookup_opacity(lookup,
                        opacity,
                        species=None,
                        vmrpert=0,
                        ax=None,
                        oneline=False,
                        total=False):
    """Plot the opacity for one or more species of an ARTS lookup table.

    Parameters:
        lookup (pyarts.classes.GasAbsLookup): ARTS lookup table.
        opacity (ndarray): Opacity per species in lookup table as generated by
            `calc_opacity_from_lookup`.
        species (pyarts.classes.ArrayOfArrayOfSpeciesTag), optional):
            ARTS species tags. If none is given, plots all species in the lookup
            table for the given vmr perturbation.
        vmrpert (int): Index of vmr perturbation for nonlinear species to plot.
        ax (AxesSubplot, optional): Axes to plot in.
        oneline (bool, optional): Draw a line where opacity == 1.
        total (bool, optional): Additionally plot the sum of opacities of all
            species.
    """
    if ax is None:
        ax = plt.gca()

    ax.set_yscale('log')
    if species is None:
        species = lookup.specs

    for tag in species:
        ax.plot(lookup.f_grid.data,
                opacity[_get_lookup_species_index(lookup, tag, vmrpert), :],
                label=',\n'.join([str(t.as_string) for t in tag]))
    if oneline:
        ax.plot(lookup.f_grid.data,
                np.ones_like(lookup.f_grid.data),
                linewidth=1,
                linestyle='--',
                color='k')
    if total:
        if lookup.nonlinspecs is not None:
            speciescount = _calc_lookup_species_count(lookup)
            spindex = np.cumsum(speciescount)
            spindex[1:] = spindex[0:-1]
            spindex[0] = 0
            spindex[[int(s) for s in lookup.nonlinspecs]] += vmrpert
            o = opacity[spindex]
        else:
            o = opacity
        ax.plot(lookup.f_grid, np.sum(o, axis=0), linewidth=1, color='k')

    if len(species) > 1:
        ax.legend(fontsize='xx-small', frameon=False)
    else:
        ax.set_title(',\n'.join(
            re.sub(r'(-\*)+$', '', str(s.as_string)) for s in species[0]),
                     y=1. - len(species[0]) * 0.05,
                     fontsize='xx-small')

    def formatter(x, pos):
        return rf'${x/1e9:g}$'

    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_minor_formatter(formatter)
    ax.tick_params(axis='both', which='major', labelsize='xx-small')
    ax.tick_params(axis='both', which='minor', labelsize='xx-small')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    return ax


def calc_opacity_from_lookup(lookup,
                             z=None,
                             g=constants.g,
                             r=_gas_constant_dry_air,
                             tpert=0):
    """Calculate the opacity from an ARTS lookup table.

    Parameters:
        lookup (pyarts.classes.GasAbsLookup): ARTS lookup table.
        z (ndarray, optional): Altitude profile. If not given, the layer
            thicknesses are calculated based on the hypsometric formula.
        g (float, optional): Gravity constant. Default uses Earth's gravity.
        r (float, optional): Gas constant for dry air. Default uses constant
            for Earth.
        tpert (int, optional): Index of temperature perturbation to plot.

    Returns:
        ndarray: Opacity per species in lookup table.
    """
    speciescount = _calc_lookup_species_count(lookup)
    vmrs = (np.repeat(lookup.vmrs.data, speciescount, axis=0)
            if lookup.nonlinspecs is not None else lookup.vmrs.data)

    ni = (lookup.p_grid.data * vmrs / lookup.t_ref.data /
          constants.Boltzmann).reshape(np.sum(speciescount), 1,
                                       len(lookup.p_grid.data))

    alpha = ni * lookup.xsec.data[tpert, :, :, :]

    if z is not None:
        z = interp1d(z.grids[0], z.data[:, 0, 0])(lookup.p_grid.data)
    else:
        # Calculate z from hypsometric formula
        pgrid = lookup.p_grid.data
        z = [
            r * t / g * np.log(p1 / p2)
            for p1, p2, t in zip(pgrid[:-1], pgrid[1:], (
                lookup.t_ref.data[1:] + lookup.t_ref.data[:-1]) / 2.)
        ]
        z = np.cumsum(z)
        p = (pgrid[1:] + pgrid[:-1]) / 2.
        z = interp1d(p, z, fill_value='extrapolate')(lookup.p_grid.data)

    return np.vstack([np.trapz(ialpha, z, axis=1) for ialpha in alpha])


def plot_arts_lookup(lookup,
                     opacity=True,
                     z=None,
                     g=constants.g,
                     r=_gas_constant_dry_air,
                     tpert=0,
                     vmrpert=0,
                     pressures=None,
                     cols=3,
                     species=None):
    """Visualize an ARTS lookup table.

    Plots the opacity or the absorption cross sections based on an
    ARTS lookup table.

    Parameters:
        lookup (pyarts.classes.GasAbsLookup): ARTS lookup table object.
        opacity (bool, optional): Set to False to plot the absorption cross
            sections.
        z (ndarray, optional): Altitude profile. Optional input for opacity
            calculation. If not given, the layer thicknesses are calculated
            based on the hypsometric formula.
        g (float, optional): Gravity constant. Uses Earth's gravity by default.
        r (float, optional): Gas constant for dry air.
            Uses constant for Earth by default.
        tpert (int, optional): Index of temperature perturbation to plot.
        vmrpert (int, optional): Index of vmr perturbation for nonlinear
            species to plot.
        pressures (ndarray(int), optional): Pressure levels to plot. If not
            given, up to 6 pressure levels are selected.
        cols (int, optional): Species to plot per row.

    Returns:
        matplotlib.figure.Figure, ndarray(AxesSubplot):
            Matplotlib Figure and Axes objects.

    Examples:

    .. plot::
        :include-source:

        from os.path import join, dirname
        import matplotlib.pyplot as plt
        import typhon as ty
        import pyarts

        lookup_file = join(dirname(ty.__file__), 'tests', 'plots',
                           'reference', 'abs_lookup_small.xml')
        fig, ax = ty.plots.plot_arts_lookup(pyarts.xml.load(lookup_file))

        fig.suptitle('Lookup table opacities')
        fig.subplots_adjust(top=0.88)
        plt.show()

    .. plot::
        :include-source:

        from os.path import join, dirname
        import matplotlib.pyplot as plt
        import typhon as ty
        import pyarts

        lookup_file = join(dirname(ty.__file__), 'tests', 'plots',
                           'reference', 'abs_lookup_small.xml')
        fig, ax = ty.plots.plot_arts_lookup(pyarts.xml.load(lookup_file),
                                            opacity=False)

        fig.suptitle('Lookup table absorption cross sections [m$^2$]')
        fig.subplots_adjust(top=0.88)
        plt.show()

    """
    if species is None:
        species = lookup.specs

    rows, cols, fig, ax = _setup_lookup_figure(lookup, cols, species)

    if opacity:
        lookup_opacity = calc_opacity_from_lookup(lookup, z, g, r, tpert)
    
    for cax, spec in zip_longest(
            ax.flatten() if len(ax.shape) == 2 else ax.reshape(ax.size, 1),
            species.data):
        if spec is None:
            cax.axis('off')
            continue

        if opacity:
            plot_lookup_opacity(lookup,
                                lookup_opacity,
                                vmrpert=vmrpert,
                                oneline=True,
                                total=True,
                                species=[spec],
                                ax=cax)
        else:
            psize = lookup.p_grid.data.size
            if pressures is not None:
                ipressures = [
                    np.abs(lookup.p_grid.data - p).argmin() for p in pressures
                ]
            else:
                ipressures = (lookup.p_grid.data.size - 1 -
                              (range(psize) if psize <= 5 else np.linspace(
                                  0,
                                  lookup.p_grid.data.size,
                                  num=6,
                                  endpoint=False,
                                  dtype=int)))
            plot_lookup_xsec(lookup,
                             ipressures,
                             species=[spec],
                             ax=cax,
                             tpert=tpert,
                             vmrpert=vmrpert)

    if opacity:
        _add_opacity_legend(ax[-1, 0])
    else:
        _add_xsec_legend(lookup, ipressures, ax[-1, 0])

    for cax in ax[-2, :]:
        cax.set_xlabel('Frequency [GHz]', fontsize='xx-small')

    return fig, ax
