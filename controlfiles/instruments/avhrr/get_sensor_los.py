#!/usr/bin/env python

"""Generates avhrr.sensor_los.xml
"""

from numpy import linspace, newaxis

import PyARTS.artsXML

sensor_los = 180 + linspace(-55.4, 55.4, 2048)[:1024]
PyARTS.artsXML.save(sensor_los[:, newaxis], "avhrr.sensor_los.xml")
