"""
This file will showcase how you can load line data from arts-cat-data
into the workspace and generate a plot of the line strengths.
"""

import os
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pyarts3 as pa
from pyarts3.arts.convert import freq2kaycm

# Download catalogs
pa.data.download()

# Initialize ARTS
ws = pa.Workspace()

# Isotopologue
isotopologue = "CO2-626"
# Frequency range
f_min = 18.7e12
f_max = 21.5e12
# Reference temperature for line strength
T0 = 296

# Load line catalog
ws.absorption_speciesSet(species=[isotopologue])
ws.absorption_bandsReadSpeciesSplitCatalog(basename="lines/")

# Throw away lines outside the frequency range
ws.absorption_bandsSelectFrequencyByLine(fmin=f_min, fmax=f_max)

# Bail out if no lines were found inside the frequency range
if len(ws.absorption_bands) == 0:
    raise RuntimeError("No lines found in frequency range")

# Get the intensity of the strongest line in each band
max_strengths = [
    np.max([line.hitran_s(isotopologue, T0=T0) for line in b.lines])
    for b in ws.absorption_bands.values()
]

# Sort bands by descending maximum line strength
indices = np.argsort(max_strengths)[::-1]
keys = list(ws.absorption_bands.keys())
sorted_bands = {keys[i]: ws.absorption_bands[keys[i]] for i in indices}

# Plot
fig, ax = plt.subplots()

# Plot the two bands with the highest line strengths
for i, (qi, band) in enumerate(list(sorted_bands.items())[:2]):
    # Line center frequencies
    pos = np.array([line.f0 for line in band.lines])
    strengths = [line.hitran_s(isotopologue, T0=T0) for line in band.lines]

    # Print band and line information
    band_name = f"Band {i + 1}"
    print(f"{band_name} info:")
    print(" ", "\n  ".join(textwrap.wrap(f"{qi:h}", width=80)))
    print(f"{band_name} lines:")
    for line, strength in zip(band.lines, strengths):
        line_info = f"{line:h}".split(";")
        strength = f"Line strength: {strength:.4g}"
        print(f"  {line_info[0].ljust(25)} {strength.ljust(25)} {line_info[-2]}")
    print()

    color = ax._get_lines.get_next_color()
    ax.vlines(
        pos / 1e12,
        0,
        strengths,
        color=color,
        label="\n".join(textwrap.wrap(f"{band_name}: {qi:h}", width=40)),
    )

ax.legend(fontsize="x-small")
ax.set_xlim(left=f_min / 1e12, right=f_max / 1e12)
ax.set_xlabel("Frequency [THz]")
ax_wn = ax.twiny()
ax_wn.set_xlim(freq2kaycm(f_min), freq2kaycm(f_max))
ax_wn.set_xlabel("Wavenumber [cm$^{-1}$]")
ax.set_ylabel("Hitran line strength")
ax.set_title(isotopologue)
ax.set_ylim(bottom=0)

if "ARTS_HEADLESS" not in os.environ:
    plt.show()
