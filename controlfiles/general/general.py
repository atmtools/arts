# DEFINITIONS:  -*-sh-*-
#
############
# General ARTS defaults
#
# The basic philosophy here is, that general.arts shall always be included by an
# ARTS controlfile. general.arts contains such settings, which are (very often)
# necessary, but on the other hand the user rarely wants/needs to reset (or is
# discouraged to do so).
# The settings done below can be divided into 3 categories:
# 1) settings that serve (mainly) as initialisations
# 2) settings the user is highly discouraged to change
# 3) settings that define some meaningful, widely valid/applicable defaults
#
############
#
# Authors: Stefan Buehler, Patrick Eriksson, Oliver Lemke, Jana Mendrok
#

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
############
# PART 1 - following settings serve (mainly) as initialisations
############
#
# Initialize verbosity levels
#
ws.verbosityInit()
#
# Semi-mandatory variables associated with scattering calculations
#
ws.ArrayOfStringSet(ws.scat_species, [])
ws.MatrixSet(ws.particle_masses, array([], shape=(1, 0), dtype=float64))
ws.Tensor4SetConstant(ws.particle_bulkprop_field, 0, 0, 0, 0, 0.0)
ws.ArrayOfStringSet(ws.particle_bulkprop_names, [])
ws.Touch(ws.dpnd_field_dx)
#
# Semi-mandatory variables associated with surface_props_data
#
ws.Tensor3SetConstant(ws.surface_props_data, 0, 0, 0, 0.0)
ws.ArrayOfStringSet(ws.surface_props_names, [])
#
# Default is that no transmitter is involved
#
ws.MatrixSet(ws.transmitter_pos, array([], shape=(1, 0), dtype=float64))
#
# Default assmption is that the sensor cause no Doppler effect
#
ws.NumericSet(ws.rte_alonglos_v, 0.0)
#
# No auxiliary variables as default
#
ws.ArrayOfStringSet(ws.iy_aux_vars, [])
#
# Wind and magnetic fields
# (all components set to be empty, shorthand for no winds/magnetic/pnd field)
#
ws.Tensor3SetConstant(ws.wind_u_field, 0, 0, 0, 0.0)
ws.Tensor3SetConstant(ws.wind_v_field, 0, 0, 0, 0.0)
ws.Tensor3SetConstant(ws.wind_w_field, 0, 0, 0, 0.0)
ws.Tensor3SetConstant(ws.mag_u_field, 0, 0, 0, 0.0)
ws.Tensor3SetConstant(ws.mag_v_field, 0, 0, 0, 0.0)
ws.Tensor3SetConstant(ws.mag_w_field, 0, 0, 0, 0.0)
############
# PART 2 - following settings are highly recommended to be NOT changed
############
#
# Set default interpolation orders for absorption lookup table. Do not
# mess with these values, unless you know what you are doing!
#
ws.IndexSet(ws.abs_p_interp_order, 5)
ws.IndexSet(ws.abs_t_interp_order, 7)
ws.IndexSet(ws.abs_nls_interp_order, 5)
# abs_f_interp_order should normally be set to 0. However, if you are doing
# calculations with Doppler shift (e.g., wind, planet rotation, or satellite
# movement considered in forward calculation), you NEED to set it to >=l. This
# regardless of whether you use absorption lookup tables or on-the-fly
# calculation of absorption. In case of on-the-fly it will have no practical
# effects. In case of lookup tables, you should choose abs_f_interp_order
# depending on the frequency interpolation scheme you want.
ws.IndexSet(ws.abs_f_interp_order, 0)
#
# Variable for calculation of propagation paths:
#
# You should not change the value of ppath_inside_cloudbox_do, if you don't
# know exactly what you are doing!)
ws.FlagOff(ws.ppath_inside_cloudbox_do)
############
# PART 3 - following settings define some meaningful defaults
############
#
# Default output format
#
ws.output_file_formatSetAscii()
#
# No unit conversion
#
ws.StringSet(ws.iy_unit, "1")
#
# Batch calculations start at index 0 by default
#
ws.IndexSet(ws.ybatch_start, 0)
#
# Variables for calculation of propagation paths:
#
# The value for ppath_lmax of 10e3 is OK for limb sounding, and
# also for down-looking if not very accurate results are demanded.
#
ws.NumericSet(ws.ppath_lmax, 10000.0)
#
# The value for ppath_lraytrace of 1e3 should be OK in general for passive
# observation, but lower values are needed for simulating radio links
#
ws.NumericSet(ws.ppath_lraytrace, 1000.0)
#
# Various geo-positioning
#
# Default is to leave the geo-pos undefined
ws.VectorSet(ws.lat_true, array([], dtype=float64))
ws.VectorSet(ws.lon_true, array([], dtype=float64))


@arts_agenda
def geo_pos_agenda(ws):
    ws.Ignore(ws.ppath)
    ws.VectorSet(ws.geo_pos, array([], dtype=float64))


ws.geo_pos_agenda = geo_pos_agenda

#
# iy_id = 0 means non-identified calculation
#
ws.IndexSet(ws.iy_id, 0)
#
# MC
#
ws.IndexSet(ws.mc_min_iter, 100)
ws.NumericSet(ws.mc_taustep_limit, 0.1)
#
# Turn off non-LTE calculations by default
#
ws.nlteOff()
#
# Use built-in Partition functions per default
#
ws.partition_functionsInitFromBuiltin()
#
# Only one option so far for water_p_eq_agenda
#
@arts_agenda
def water_p_eq_agenda(ws):
    ws.water_p_eq_fieldMK05()


ws.water_p_eq_agenda = water_p_eq_agenda

ws.IndexSet(ws.scat_data_checked, 0)
