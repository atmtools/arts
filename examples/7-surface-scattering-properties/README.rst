These examples work with the surface scattering framework: a map of scattering
models keyed by surface property tags, used as the boundary condition of the
*spectral_rad_surface_agenda*.

The surface field carries one mask per tag, and the agenda resolves the model of
a hit point from the mask that is set there.  Two ways of using that are shown
here:

* ``surface-scattering-framework.py`` only builds the objects of the framework,
  the reflectivity fields, the scatterers, and the map holding them, without
  running any radiative transfer;
* ``1-surface-scattering-checkerboard.py`` uses them as the lower boundary of a
  full simulation and plots what a sensor in a high orbit sees of a surface
  painted in two spectrally opposite colours;
* ``2-surface-scattering-checkerboard-atm-lut.py`` repeats the same scene
  through a tropical clearsky atmosphere, with the absorption taken from a
  lookup table and the sky lit up by first-order solar Rayleigh scattering;
* ``3-surface-scattering-checkerboard-atm-lut-ir.py`` moves the same scene to
  the thermal infrared -- no sun, two window wavenumbers, emissivity from
  Kirchhoff's law, and a latitudinal surface temperature pattern.
