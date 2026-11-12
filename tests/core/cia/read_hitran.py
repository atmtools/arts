"""Reading HITRAN CIA files with abs_cia_dataReadFromCIA

The method looks for the files of a CIA pair in the folders Main-Folder/<pair>
and Alternate-Folder/<pair> of a catalog.  A file of two temperatures, each
with three wavenumbers, is written there in the HITRAN format and read back in
SI units.
"""

import os
import tempfile

import numpy as np
import pyarts3 as pyarts

WAVENUMBERS = np.array([10.0, 20.0, 30.0])  # [cm-1]
TEMPERATURES = np.array([100.0, 200.0])  # [K]
CROSS_SECTIONS = np.array([[1.0e-46, 2.0e-46, 3.0e-46], [4.0e-46, 5.0e-46, 6.0e-46]])  # [cm5 molecule-2] by temperature


def header(temperature):
    """The HITRAN CIA header of a data set: A20, F10.3, F10.3, I7, F7.1, E10.3, F6.3, A27, I3"""
    return (
        f"{'H2-He':>20s}{WAVENUMBERS[0]:10.3f}{WAVENUMBERS[-1]:10.3f}{len(WAVENUMBERS):7d}{temperature:7.1f}"
        f"{CROSS_SECTIONS.max():10.3e}{0.0:6.3f}{'':27s}{1:3d}"
    )


with tempfile.TemporaryDirectory() as catalog:
    folder = os.path.join(catalog, "Alternate-Folder", "H2-He")
    os.makedirs(folder)
    with open(os.path.join(folder, "H2-He_test.cia"), "w") as file:
        for temperature, cross_sections in zip(TEMPERATURES, CROSS_SECTIONS):
            file.write(header(temperature) + "\n")
            for wavenumber, cross_section in zip(WAVENUMBERS, cross_sections):
                file.write(f"{wavenumber:10.4f} {cross_section:10.3e}\n")

    ws = pyarts.Workspace()
    ws.abs_speciesSet(species=["H2-CIA-He"])
    ws.abs_cia_dataReadFromCIA(catalogpath=catalog)

(record,) = ws.abs_cia_data.values()
(field,) = record.data
assert np.allclose(field.grids[0], WAVENUMBERS * 100 * pyarts.arts.constants.c, rtol=1e-15, atol=0)
assert np.array_equal(field.grids[1], TEMPERATURES)
assert np.allclose(field.data, CROSS_SECTIONS.T * 1e-10, rtol=1e-15, atol=0)
