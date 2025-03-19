.. _Sec Sensor description:

Sensor description
##################

The sensor modelling in ARTS is meant to simply describe how we take 
Stokes vectors from one or several spectral raditive transfer simulations
and turn them into scalar measurements.

The idea is that each sensor has a set of positions and line-of-sights
where it takes measurements (its observation geometry), a set of frequencies
that it samples (its backend sampling),
and a sensitivity to some atmospheric state of polarization.

The main variable that describes 
the sensor setup is :attr:`~pyarts.workspace.Workspace.measurement_sensor`.
It is a list of sensor observation elements.
Every sensor observation element (:class:`~pyarts.arts.SensorObsel`) has a
shared frequency grid and shared geometry.

Sensor observation elements
***************************

A sensor observation element describes how the sensor samples a single
scalar measurement. 

The equation for a single scalar measurement is:

.. math::

    y =
    \sum_g \sum_f \left(\vec{w}_{g, f} \cdot \vec{I}_{g, f}\right) + \epsilon =
    y_0 + \epsilon,

where :math:`y` is the scalar measurement, :math:`\vec{w}_{g, f}` is the polarized weight,
:math:`\vec{I}_{g, f}` is the weight-corresponding Stokes vector, and
:math:`\epsilon` is the systematic error that is beyond forward simulations.
The :math:`y_0`-term
is introduced as a way to separate the forward model from the error estimate.
The indices :math:`g` and :math:`f` are for the geometry
and frequency, respectively.

.. tip::

    It is generally a good idea if the sum of the weights multiplied by the
    Stokes vector describing your radiation sensitivity is close to 1.

Backend sampling
----------------

The frequency grid determines what frequencies are sampled by the sensor.
The frequency grid is shared between all sensor observation elements belonging
to the same sensor.  The weight 

Antenna pattern
---------------

Antenna patterns describe the sensitivity of the sensor to the incoming
radiation.  The antenna pattern is a function of the angles between the
incoming radiation and the sensor line-of-sight.  An antenna pattern is
modelled in ARTS using pure changes in the observation geometry coupled
to a weight redistribution in the sensor observation element.

.. note::
 
    The geometry may be shared between multiple sensor observation
    elements but the weights are not.  This is useful for speeding up calculations
    when the sensor observation elements share the same frequency grid and
    are considered to be the same sensor.

Sensor systematic error
=======================

The systematic error :math:`\epsilon` is for a single sensor. It is
a free-form function.  Using builtin functions for this helps when the shape
of :math:`\epsilon` should be retrieved.  Multiple :math:`\epsilon` for a
single sensor may be added.
