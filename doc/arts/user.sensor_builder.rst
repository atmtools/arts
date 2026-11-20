.. _Sec Sensor builder:

Building a sensor
#################

A :class:`~pyarts3.arts.SensorBuilder` describes an instrument: its channels and
its antenna pattern.  Placing it at sensor positions and bore lines of sight
gives a :class:`~pyarts3.arts.Sensor`, as described in
:ref:`Sec Sensor description`, to put in
:attr:`~pyarts3.workspace.Workspace.measurement_sensors`.  The instrument is kept
apart from where it looks, so the same builder serves any geometry, and it can
be saved to and loaded from an ARTS XML file.

``examples/6-satellite-sensors/amsu-a.py`` loads an AMSU-A-like instrument from
a file and places it in a cross-track scan.

What a builder holds
********************

A builder holds:

- ``channels``: a list of :class:`~pyarts3.arts.sensors.Channel`, one
  observation element per channel, pixel of the antenna and placement; and
- ``antenna``: one antenna pattern, shared by all channels.

Both can be set after construction.  A builder can be constructed from

- a list of channels and an antenna pattern;
- a :class:`~pyarts3.arts.sensors.Spectrometer` and an antenna pattern, which
  keeps the spectrometer's channels on one common frequency grid; or
- an antenna pattern, a spectrometer and a heterodyne receiver
  (:class:`~pyarts3.arts.sensors.HeterodyneFrequencyRange`), which takes the
  spectrometer's channels through the receiver to the real frequencies.

Channels
********

A :class:`~pyarts3.arts.sensors.Channel` is a frequency response and a
polarization response.

The frequency response is a sum of pieces.  Each piece is linear between its
frequencies and zero outside them, and a piece with a single frequency is a
delta function.  The frequencies of a piece are also where the radiance is
sampled, so the sampling density is chosen when the channel is defined.  The
channel measures the integral of the response times the radiance over the
integral of the response; :meth:`~pyarts3.arts.sensors.Channel.quadrature`
gives the weights.  Frequencies are in Hz.

Channels can be made from response pieces directly or with these helpers:

- :class:`~pyarts3.arts.sensors.BoxChannel`: a uniform response over an
  interval;
- :class:`~pyarts3.arts.sensors.DiracChannel`: a single frequency; and
- :class:`~pyarts3.arts.sensors.GaussianChannel`: a Gaussian response.

Their frequencies are used as they are, so they can be real frequencies, or
frequencies on the intermediate-frequency axis of a receiver.

A :class:`~pyarts3.arts.sensors.Spectrometer` is a set of channels on shared
frequencies: one channel shifted to a list of offsets, Dirac channels at a list
of offsets, or a given list of channels.

A :class:`~pyarts3.arts.sensors.HeterodyneFrequencyRange` describes a receiver
chain as a sequence of filters and mixer stages.  Its
:meth:`~pyarts3.arts.sensors.HeterodyneFrequencyRange.channel_response` turns a
channel defined on the final intermediate-frequency axis into the channel seen
at the real frequencies.  Each response piece becomes one piece per sideband
that passes it, so the sidebands are summed.  The polarization is kept.

Polarization
************

A :class:`~pyarts3.arts.sensors.Polarization` has a state and a scan rotation.

The state is the Stokes vector of weights at the boresight, in the polarization
basis of the received radiation, that is of its propagation direction, which is
opposite to the line of sight (see :ref:`Sec Polarization basis`).  It can be
given as a Stokes vector or as the name of a
:class:`~pyarts3.arts.PolarizationChoice`, so ``"Iv"`` receives vertical and
``"Ih"`` horizontal polarization.

The scan rotation turns the linear part of the state about the ray as the
instrument scans, by that many degrees per degree of scan angle.  It is 0 when
the polarization does not depend on the pointing.  A cross-track scanner with a
rotating reflector has a scan rotation of +1 or -1, and its state is then the
polarization at nadir.  The sign belongs to the instrument and only matters when
the observed U is not zero.  :meth:`~pyarts3.arts.sensors.Polarization.at` gives
the state at a scan angle.  A channel with a nonzero scan rotation can only be
placed with a scan geometry, see `Placing the instrument`_.

The antenna pattern carries the boresight polarization to each of its lines of
sight along the great circle from the boresight, and expresses it in the
polarization basis of that line of sight.  The patterns have no
cross-polarization.

Antenna patterns
****************

``antenna`` holds one antenna pattern.  An antenna pattern describes the layout
of its pixels and the beam of each pixel, and gives one observation element per
pixel for each channel.  The layout, its ``pixels``, is a
:class:`~pyarts3.arts.sensors.PixelGrid`: rows and columns of angular offsets
from the boresight, with the observation elements ordered by row and then by
column.  A camera has a grid of pixels, and the other kinds are a single beam,
one pixel at the boresight.

Each kind of pattern is its own class in ``pyarts3.arts.sensors``, and
``antenna`` accepts any of them.  The
documentation of :attr:`~pyarts3.arts.SensorBuilder.antenna` lists the kinds
that exist.  More kinds can be added, so code that depends on the kind should
check the type of the pattern rather than assume one.  A saved builder records
the kind of its pattern.

The kinds are:

- :class:`~pyarts3.arts.sensors.PencilBeamAntenna`: only the bore line of sight.
- :class:`~pyarts3.arts.sensors.GaussianAntenna`: a gain
  :math:`\exp\left(-\frac{1}{2}(\theta / \sigma)^2\right)` of the angular
  offset :math:`\theta` from the boresight.  It holds the standard deviation
  :math:`\sigma`, the zenith offsets and the azimuth offsets at which the gain
  is sampled, all in degrees.  The gain does not depend on frequency.
- :class:`~pyarts3.arts.sensors.GaussianAiryAntenna`: a Gaussian gain whose width
  is proportional to the wavelength over the aperture diameter, so it narrows
  with frequency.  It holds the aperture diameter in m and the zenith and
  azimuth offsets.  The channel frequencies must be positive.
- :class:`~pyarts3.arts.sensors.GriddedAntennaPattern`: any gain, given as an
  :class:`~pyarts3.arts.sensors.AntennaPatternField` on zenith and azimuth
  offsets from the boresight.
- :class:`~pyarts3.arts.sensors.CameraAntenna`: a grid of pixels around the
  optical axis, the boresight, each with the same beam, which is one of the
  single-beam kinds above.  Its pixel grid is in CCD coordinates: rows towards
  increasing zenith and, because the lens inverts the image, columns towards
  decreasing azimuth.  A thin-lens camera is made from its CCD, focal length
  and focus distance with
  :meth:`~pyarts3.arts.sensors.CameraAntenna.thin_lens`.  Each pixel receives
  the boresight polarization carried to it along the great circle.

For all kinds that sample more than the boresight, each sample stands for a
cell of the sphere around the boresight: the ring between the midpoints to its
neighbouring zenith offsets, shared between its azimuths.  The outermost ring
reaches as far beyond its offset as half the distance to its inner neighbour.
A zenith offset of 0
is a single sample, the cap inside the first midpoint, with the mean gain over
its azimuths.  The weight of a line of sight is its gain times the solid angle
of its cell, normalized to sum to one at each frequency.

Each pattern can also be called with a channel and a geometry, which gives a
sensor with one observation element per pixel, for inspection.

Placing the instrument
**********************

Calling the builder with positions, bore lines of sight and the reference
ellipsoid pairs each position with its line of sight, one placement of the
instrument each.  It returns one :class:`~pyarts3.arts.Sensor`: the frequencies
of all channels, the lines of sight that the antenna samples at every
placement, and the observation elements ordered by placement, then by pixel of
the antenna, then by channel.  Its meta has these axes: placement by channel
for an antenna with one pixel, and placement by rows by columns of pixels by
channel for a camera.  The sensor is named by the optional ``name``, or after
its number of channels without one.
Positions and lines of sight are as in :class:`~pyarts3.arts.SensorPosLos`.

:meth:`~pyarts3.arts.SensorBuilder.cross_track` places a cross-track scanner.
Each position has a heading, its flight direction in degrees clockwise from
north, and scans the same signed scan angles across that direction.  A scan
angle :math:`s` is the nadir angle of the boresight, positive to the right of
the flight direction, so the bore line of sight is
:math:`[180 - |s|, \mathrm{heading} + 90]` for :math:`s \geq 0` and
:math:`[180 - |s|, \mathrm{heading} - 90]` for :math:`s < 0`.  Each channel's
polarization is taken at :math:`s`.  Each position and scan angle is one
placement, ordered by position and then scan angle.  The antenna must have one
pixel.

All the elements of the sensor share its frequency grid and its lines of
sight, so each line of sight is simulated once over the frequencies of all
channels, see :ref:`Sec Sensor description`.

Saving and loading
******************

A builder is saved with ``savexml`` and loaded with
:meth:`~pyarts3.arts.SensorBuilder.fromxml` or :func:`pyarts3.xml.load`.  The
file holds the channels with their response pieces and polarization, and the
antenna pattern with its kind and parameters.  It does not hold any geometry.
