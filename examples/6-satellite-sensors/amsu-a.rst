AMSU-A, the Advanced Microwave Sounding Unit-A, is a cross-track scanning
microwave radiometer on NOAA-15 to -19, Aqua and MetOp-A to -C.  Channels 3
to 14 sound temperature in the 50-58 GHz oxygen band, and channels 1, 2 and 15
see the surface, water vapour and clouds at 23.8, 31.4 and 89 GHz.

``amsu-a.xml`` holds a :class:`~pyarts3.arts.SensorBuilder` with

* the pass bands of each channel, idealised as boxcars around the nominal
  frequencies and local-oscillator offsets,
* the polarization of each channel: vertical at nadir for channels 1-4, 7 and
  15, horizontal for the others, turning with the scan angle as the reflector
  rotates,
* a pencil-beam antenna.

The scan geometry is not part of the sensor.  It is given when the sensor is
placed, here as one position scanning at 50 degrees.  The file does not hold
the real instrument's antenna pattern, noise or calibration.
