#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 10:44:35 2026

@author: u242031
"""
import numpy as np
import pyarts3.arts as paa3


Lambertian=paa3.LambertianSurfaceScatterer()

surface_models=paa3.MapOfSurfaceScatteringModel()

surface_models.add('land',Lambertian)