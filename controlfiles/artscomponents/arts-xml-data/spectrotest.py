import numpy as N
import pylab as p
import os

# import pdb
from matplotlib.gridspec import GridSpec as gc


def get_data(fa3h=None, fa4h=None, fa4p=None):
    """
  Reads data (cross sections, frequency grid, pressure grid, absorption species)
   from ARTS lookup table files.

  Parameters
  ----------
  fa3h: string
      name of ARTS3 catalogue file from HITRAN data
  fa4h: string
      name of ARTS4 catalogue file from HITRAN data
  fa4p:
      name of ARTS4 catalogue file from toolbox spectroscopic line data
  NOTE 1: not all these files have to be given. 'None' will create empty
   cross-section output. If none of them is given, also the other output will be
   empty.
  NOTE 2: all of them might be used for other ARTS LUT data, too. Basically we
   can process 3 LUT files, whatever they contain (but all of them should be for
   identical frequency, pressure, and species grids!!!).

  Returns
  -------
  freq: array of float (dim: nfreq)
      Frequency grid
  pres: array of float (dim: npres)
      Pressure grid
  absspec: list of string (dim: nspec)
      Species name tags
  axs_a3h: array of float (dim: (ntemp,nspec,nfreq,npres))
      absorption cross section tensor from ARTS3 HITRAN data
  axs_a4h: array of float (dim: (ntemp,nspec,nfreq,npres))
      absorption cross section tensor from ARTS4 HITRAN data
  axs_a4p: array of float (dim: (ntemp,nspec,nfreq,npres))
      absorption cross section tensor from ARTS4 Toolbox data
  NOTE: ntemp is assumed to be 1.
  """
    # pdb.set_trace()
    if fa3h is not None:
        f = open(fa3h)
    elif fa4h is not None:
        f = open(fa4h)
    elif fa4p is not None:
        f = open(fa4p)
    else:
        absspec = []
        freq = N.array([])
        pres = N.copy(freq)
        axs_a3h = N.copy(freq)
        axs_a4h = N.copy(freq)
        axs_a4p = N.copy(freq)
        return absspec, freq, pres, axs_a3h, axs_a4h, axs_a4p

    lines = f.readlines()
    f.close()

    # this works only with one tag per species. and it neglects possible isotopologue subtags.
    # absspec=[x for x in lines if '<SpeciesTag>' in x]
    # for i,l in enumerate(absspec): absspec[i]=l.partition('"')[2].partition('-')[0]

    # so, instead we're trying to be smarter (looks less beautiful but doing what i want):
    i1 = first_substring(lines, "ArrayOfSpeciesTag")
    nabsspec = N.int(lines[i1].rpartition('nelem="')[2].rpartition('"')[0])
    absspec = []
    for i in N.arange(nabsspec):
        i1 += 1
        i2 = first_substring(lines[i1:], 'type="SpeciesTag"')
        assert i2 >= 0, (
            'Looking for %i "ArrayOfSpeciesTag" entries, but could find only %s.'
            + " Aborting.\n" % (nabsspec, i)
        )
        i1 = i1 + i2
        ntags = N.int(lines[i1].rpartition('nelem="')[2].rpartition('"')[0])
        isostrings = ""
        for j in N.arange(ntags) + 1:
            isostring = lines[i1 + j].partition('"')[2]
            if "*" in isostring:
                isostring = isostring.partition("-*")[0]
            else:
                isostring = isostring.partition('"')[0]
            isostrings = isostrings + "," + isostring
        # removing the comma we added for loop simplicity before the first tag
        absspec.append(isostrings[1:])
    print(absspec)

    ind = first_substring(lines, "FrequencyGrid")
    nfreq = N.int(lines[ind].partition('nelem="')[-1].partition('"')[0])
    freq = N.zeros(nfreq)
    for i in N.arange(nfreq):
        freq[i] = N.float(lines[ind + i + 1])

    ind = first_substring(lines, "PressureGrid")
    npres = N.int(lines[ind].partition('nelem="')[-1].partition('"')[0])
    pres = N.zeros(npres)
    for i in N.arange(npres):
        pres[i] = N.float(lines[ind + i + 1])

    ind = first_substring(lines, "AbsorptionCrossSections")
    if fa3h is None:
        axs_a3h = N.array([])
    else:
        axs_a3h = N.loadtxt(fa3h, skiprows=ind + 1, comments="<").reshape(
            (-1, nabsspec, nfreq, npres)
        )
    if fa4h is None:
        axs_a4h = N.array([])
    else:
        axs_a4h = N.loadtxt(fa4h, skiprows=ind + 1, comments="<").reshape(
            (-1, nabsspec, nfreq, npres)
        )
    if fa4p is None:
        axs_a4p = N.array([])
    else:
        axs_a4p = N.loadtxt(fa4p, skiprows=ind + 1, comments="<").reshape(
            (-1, nabsspec, nfreq, npres)
        )

    return absspec, freq, pres, axs_a3h, axs_a4h, axs_a4p


def fast_stats(absspec, axs1, axs2, axs3=None, header=1):
    """
  Calculates and prints out some deviation statistics (max and mean relative
   deviation) between two (or three) sets of absorption cross sections.

  Parameters
  ----------
  absspec: string
      Species name tags corresponding to species in axs1/2 (2nd dimension)
  axs1/axs2/axs3: array of float (dim: (ntemp,nspec,nfreq,npres))
      The absorption cross section tensors. Sets axs1 and axs2 are mandatory. If
       axs3 is present, it is compared to the axs1.    

  Returns
  -------
  no returns
  """
    if axs3 is not None:
        if header:
            print(
                "%10s %15s %15s %15s %15s"
                % (
                    "Species",
                    "max(ax1-ax2)",
                    "mean(ax1-ax2)",
                    "max(ax1-ax3)",
                    "mean(ax1-ax3)",
                )
            )
        for i in N.arange(axs1.shape[1]):
            print(
                "%10s %15.1e %15.1e %15.1e %15.1e"
                % (
                    absspec[i],
                    N.nanmax(
                        2.0
                        * abs(axs1[0, i, :, :] - axs2[0, i, :, :])
                        / (axs1[0, i, :, :] + axs2[0, i, :, :])
                    ),
                    N.nansum(
                        2.0
                        * abs(axs1[0, i, :, :] - axs2[0, i, :, :])
                        / (axs1[0, i, :, :] + axs2[0, i, :, :])
                    )
                    / N.isfinite(1.0 / (axs1[0, i, :, :] + axs2[0, i, :, :])).sum(),
                    N.nanmax(
                        2.0
                        * abs(axs1[0, i, :, :] - axs3[0, i, :, :])
                        / (axs1[0, i, :, :] + axs3[0, i, :, :])
                    ),
                    N.nansum(
                        2.0
                        * abs(axs1[0, i, :, :] - axs3[0, i, :, :])
                        / (axs1[0, i, :, :] + axs3[0, i, :, :])
                    )
                    / N.isfinite(1.0 / (axs1[0, i, :, :] + axs3[0, i, :, :])).sum(),
                )
            )
    else:
        if header:
            print("%10s %15s %15s" % ("Species", "max(ax1-ax2)", "mean(ax1-ax2)"))
        for i in N.arange(axs1.shape[1]):
            print(
                "%10s %15.1e %15.1e"
                % (
                    absspec[i],
                    N.nanmax(
                        2.0
                        * abs(axs1[0, i, :, :] - axs2[0, i, :, :])
                        / (axs1[0, i, :, :] + axs2[0, i, :, :])
                    ),
                    N.nansum(
                        2.0
                        * abs(axs1[0, i, :, :] - axs2[0, i, :, :])
                        / (axs1[0, i, :, :] + axs2[0, i, :, :])
                    )
                    / N.isfinite(1.0 / (axs1[0, i, :, :] + axs2[0, i, :, :])).sum(),
                )
            )


def plotdiff(
    axs1,
    axs2,
    freq,
    absspec,
    start=0,
    ende=None,
    newfig=False,
    save=False,
    specplot=1,
    title="HITRAN vs. Toolbox for ",
    casename="H-vs-TB",
):
    """
  Plots deviations of two sets of absorption cross sections.
  
  Plotting is done over frequency for all pressure levels present into one plot.
   Separate plots for separate species.
  Simultaneously does four different type of plots:
  (A) relative deviations
    (A1) log plot of abs(axs1-axs2) / 0.5*(axs1+axs2)
    (A2) linear plot of (axs1-axs2) / 0.5*(axs1+axs2)
  (B) absolute deviations
    (B1) linear plot of (axs1-axs2)
    (B2) log plot of abs(axs1-axs2)

  Parameters
  ----------
  axs1/axs2: array of float (dim: (ntemp,nspec,nfreq,npres))
      The absorption cross section tensors.
  freq: array of float (dim: nfreq)
      Frequency grid.
  absspec: list of string (dim: nspec)
      Species name tags corresponding to axs1 and axs2 2nd dimension.
  start: integer
      Index of absspec species to start with (earlier entries are ignored).
  ende: integer
      Index of last absspec species to be plotted (later entries are ignored).
  newfig: boolean
      Flag whether to open a new figure for each plot (=True; recommended when
       plotting on screen). If False, figure content will be erased before each
       new plot (recommended when saving the figure, avoiding to open plenty of
       on-screen figure windows).
  save: boolean
      Flag whether to save figures as graphic (png) files.
  specplot: boolean (or integer)
      Flag whether to make plots of the abs-xsec spectra themselves.
  title: string
      Basic string of each figure title. Will be complemented with the absspec
       tag corresponging to the respective figure.
  casename: string
      Basic string for graphics file names. Graphics are saved in figures/ with
       full name constructed from a plot type indentifier, this casename and the
       abspsec tag of the corresponding species.

  Returns
  -------
  no returns
  """
    col = ["b", "g", "r", "c", "m", "y", "k"]
    if ende is None:
        ende = axs1.shape[1] - 1
    for i in N.arange(start, ende + 1):
        if (axs1[0, i, :, :] + axs2[0, i, :, :]).max() != 0.0:
            if newfig:
                p.figure()
            else:
                p.clf()
            p.semilogy(
                freq * 1e-9,
                0.5
                * abs(axs1[0, i, :, :] - axs2[0, i, :, :])
                / (axs1[0, i, :, :] + axs2[0, i, :, :])
                + 1e-60,
            )
            p.title("%s %s" % (title, absspec[i]))
            p.xlabel("Frequency [GHz]")
            p.ylabel("relative difference [-]")
            p.tight_layout()
            if save:
                p.savefig("figures/drel_logabs_%s_%s.png" % (casename, absspec[i]))
            if newfig:
                p.figure()
            else:
                p.clf()
            p.plot(
                freq * 1e-9,
                0.5
                * (axs1[0, i, :, :] - axs2[0, i, :, :])
                / (axs1[0, i, :, :] + axs2[0, i, :, :]),
            )
            p.title("%s %s" % (title, absspec[i]))
            p.xlabel("Frequency [GHz]")
            p.ylabel("relative difference [-]")
            p.tight_layout()
            if save:
                p.savefig("figures/drel_lin_%s_%s.png" % (casename, absspec[i]))
            if newfig:
                p.figure()
            else:
                p.clf()
            p.plot(freq * 1e-9, axs1[0, i, :, :] - axs2[0, i, :, :])
            p.xlabel("Frequency [GHz]")
            p.title("%s %s" % (title, absspec[i]))
            p.ylabel("absolute difference [m2]")
            p.tight_layout()
            if save:
                p.savefig("figures/dabs_lin_%s_%s.png" % (casename, absspec[i]))
            if newfig:
                p.figure()
            else:
                p.clf()
            p.semilogy(freq * 1e-9, abs(axs1[0, i, :, :] - axs2[0, i, :, :]) + 1e-60)
            p.title("%s %s" % (title, absspec[i]))
            p.xlabel("Frequency [GHz]")
            p.ylabel("absolute difference [m2]")
            p.tight_layout()
            if save:
                p.savefig("figures/dabs_logabs_%s_%s.png" % (casename, absspec[i]))
            if specplot:
                if newfig:
                    p.figure()
                else:
                    p.clf()
                for j in N.arange(axs2.shape[-1]):
                    p.plot(
                        freq * 1e-9,
                        axs2[0, i, :, j],
                        col[j % N.size(col)] + "-",
                        linewidth=2,
                    )
                    p.plot(
                        freq * 1e-9,
                        axs1[0, i, :, j],
                        col[j % N.size(col)] + "--",
                        linewidth=2,
                    )
                p.xlabel("Frequency [GHz]")
                p.title("%s %s" % (title, absspec[i]))
                p.ylabel("absorption cross section [m2]")
                p.tight_layout()
                if save:
                    p.savefig("figures/abs_lin_%s_%s.png" % (casename, absspec[i]))
                if newfig:
                    p.figure()
                else:
                    p.clf()
                for j in N.arange(axs2.shape[-1]):
                    p.semilogy(
                        freq * 1e-9,
                        axs2[0, i, :, j],
                        col[j % N.size(col)] + "-",
                        linewidth=2,
                    )
                    p.semilogy(
                        freq * 1e-9,
                        axs1[0, i, :, j],
                        col[j % N.size(col)] + "--",
                        linewidth=2,
                    )
                p.title("%s %s" % (title, absspec[i]))
                p.xlabel("Frequency [GHz]")
                p.ylabel("absolute cross section [m2]")
                p.tight_layout()
                if save:
                    p.savefig("figures/abs_log_%s_%s.png" % (casename, absspec[i]))
        else:
            print("No abs xs !=0 (no lines? vmr=0?) for %s" % absspec[i])
        p.show()


def plotabsxsec(
    axs1,
    axs2,
    freq,
    pres,
    absspec,
    abswhich=None,
    doboth=1,
    meandelta=0,
    facdiff=0,
    newfig=True,
    save=False,
    title="Spectroscopy HITRAN vs. Toolbox: ",
    casename="H-vs-TB",
    outdir="~/projects/MicrowavePropagationToolbox/study/Validation/basics/figures/",
):
    """
  """
    sps1, sps2 = prep2to1plot()
    col = ["b", "g", "r", "c", "m", "y", "k"]
    if abswhich is None:
        abswhich = N.arange(axs1.shape[1])
    for i in abswhich:
        if (axs1[0, i, :, :] + axs2[0, i, :, :]).max() != 0.0:
            p.figure()
            ax1 = p.subplot(sps1)
            ax2 = p.subplot(sps2, sharex=ax1)
            for j in N.arange(axs2.shape[-1]):
                if meandelta:
                    ax2.plot(
                        freq * 1e-9,
                        (axs2[0, i, :, j] - axs1[0, i, :, j])
                        / (axs2[0, i, :, j] + axs1[0, i, :, j])
                        * 2e2,
                        col[j % N.size(col)] + "-",
                        linewidth=2,
                    )
                elif facdiff:
                    #           ax2.semilogy(freq*1e-9,N.maximum(axs1[0,i,:,j]/axs2[0,i,:,j],axs2[0,i,:,j]/axs1[0,i,:,j])*1e2-1e2,
                    #                   col[j%N.size(col)]+'-',linewidth=2)
                    ax2.semilogy(
                        freq * 1e-9,
                        N.maximum(
                            axs1[0, i, :, axs2.shape[-1] - j - 1]
                            / axs2[0, i, :, axs2.shape[-1] - j - 1],
                            axs2[0, i, :, axs2.shape[-1] - j - 1]
                            / axs1[0, i, :, axs2.shape[-1] - j - 1],
                        )
                        - 1.0,
                        col[(axs2.shape[-1] - j - 1) % N.size(col)] + "-",
                        linewidth=2,
                    )
                else:
                    ax2.plot(
                        freq * 1e-9,
                        axs2[0, i, :, j] / axs1[0, i, :, j] * 1e2 - 1e2,
                        col[j % N.size(col)] + "-",
                        linewidth=2,
                    )
                ax1.semilogy(
                    freq * 1e-9,
                    axs2[0, i, :, j],
                    col[j % N.size(col)] + "-",
                    linewidth=2,
                    label="%.1ePa" % pres[j],
                )
                if doboth:
                    ax1.semilogy(
                        freq * 1e-9,
                        axs1[0, i, :, j],
                        col[j % N.size(col)] + "--",
                        linewidth=2,
                    )
            ax1.set_xlabel("Frequency [GHz]")
            ax1.set_ylabel("absolute cross section [m2]")
            if facdiff:
                ax2.set_ylabel("factor difference [-]")
            else:
                ax2.set_ylabel("relative difference [%]")
            ax1.legend(loc=0)
            p.minorticks_on()
            p.title("%s %s" % (title, absspec[i]))
            p.tight_layout()
            if save:
                p.savefig(
                    "%sAbsXsec_%s_%s_xsec-log_drel-lin.png"
                    % (os.path.expanduser(outdir), casename, absspec[i])
                )
        else:
            print("No abs xs !=0 (no lines? vmr=0?) for %s" % absspec[i])
        p.show()
    return ax1, ax2


def plotabsxsec3(
    axs1,
    axs2,
    freq,
    pres,
    absspec,
    abswhich=None,
    newfig=True,
    save=False,
    title="Spectroscopy HITRAN vs. Toolbox: ",
    casename="H-vs-TB",
    outdir="~/projects/MicrowavePropagationToolbox/study/Validation/basics/figures/",
):
    """
  """
    col = ["b", "g", "r", "c", "m", "y", "k"]
    if abswhich is None:
        abswhich = N.arange(axs1.shape[1])
    for i in abswhich:
        if (axs1[0, i, :, :] + axs2[0, i, :, :]).max() != 0.0:
            fig, ax = p.subplots(3, sharex=True)
            ax1 = ax[2]
            ax2 = ax[1]
            ax3 = ax[0]
            ax3.set_title("%s %s" % (title, absspec[i]))
            for j in N.arange(axs2.shape[-1]):
                ax3.plot(
                    freq * 1e-9,
                    axs2[0, i, :, j] / axs1[0, i, :, j] * 1e2 - 1e2,
                    col[j % N.size(col)] + "-",
                    linewidth=2,
                )
                ax2.semilogy(
                    freq * 1e-9,
                    abs(axs2[0, i, :, j] - axs1[0, i, :, j]),
                    col[j % N.size(col)] + "-",
                    linewidth=2,
                )
                ax1.semilogy(
                    freq * 1e-9,
                    axs2[0, i, :, j],
                    col[j % N.size(col)] + "-",
                    linewidth=2,
                    label="%.1ePa" % pres[j],
                )
            ax1.set_xlabel("Frequency [GHz]")
            ax1.set_ylabel("abs. x-sec. [m2]")
            ax2.set_ylabel("abs. diff. [m2]")
            ax3.set_ylabel("rel. diff. [%]")
            ax1.legend(loc=0)
            ax1.minorticks_on()
            ax2.minorticks_on()
            ax3.minorticks_on()
            fig.tight_layout()
            p.show()
            if save:
                p.savefig(
                    "%sAbsXsec_%s_%s_xsec-log_dabs-log_drel-lin.png"
                    % (os.path.expanduser(outdir), casename, absspec[i])
                )
        else:
            print("No abs xs !=0 (no lines? vmr=0?) for %s" % absspec[i])
        p.show()
    return ax1, ax2, ax3


def prep2to1plot():
    """
  prepare setup for plotting window with 2 vertically stacked panels, where
   lower one is double as high as the upper one.
  """
    gridspec = gc(3, 1)
    subplotspec1 = gridspec.new_subplotspec((1, 0), 2, 1)
    subplotspec2 = gridspec.new_subplotspec((0, 0), 1, 1)
    return subplotspec1, subplotspec2


def first_substring(strings, substring):
    """
  Derives the first occurence of substring within a list of strings.
  
  Parameters
  ----------
  strings: list of strings
      list with strings that shall be searched
  substring: string
      the string whose first occurence in strings shall be determined (can be a
       substring of the entries in strings)
      
  Returns
  -------
  index: integer
      the index of the first entry in strings that contains substring
      if substring is not found within the list of strings, -1 is returned
  """
    try:
        index = (i for i, string in enumerate(strings) if substring in string).next()
    except:
        index = -1
    return index


if __name__ == "__main__":
    fa3h = "output/ARTS3-from-HITRAN08_LUT_WithFascode.xml"
    fa4h = "output/ARTS4-from-HITRAN08_LUT_WithFascode.xml"
    fa4p = "output/Perrin_LUT_WithFascode.xml"
    fa3h2 = "output/ARTS3-from-HITRAN08_LUT_NoFascode.xml"
    fa4h2 = "output/ARTS4-from-HITRAN08_LUT_NoFascode.xml"
    fa4p2 = "output/Perrin_LUT_NoFascode.xml"

    absspec, freq, pres, axs_a3h, axs_a4h, axs_a4p = get_data(fa3h, fa4h, fa4p)
    absspec2, freq2, pres2, axs_a3h2, axs_a4h2, axs_a4p2 = get_data(fa3h2, fa4h2, fa4p2)

    fast_stats(absspec, axs_a3h, axs_a4h, axs_a4p)
    fast_stats(absspec2, axs_a3h2, axs_a4h2, axs_a4p2)
