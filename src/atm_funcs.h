/* Copyright (C) 2000 Patrick Eriksson <patrick@rss.chalmers.se>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */


//==========================================================================
//=== Physical functions.
//==========================================================================

/** Calculates a blackbody radiation (the Planck function) matrix.
    Each row of the returned matrix corresponds to a frequency, while each
    column corresponds to a temperature.

    @param B       output: the blackbody radiation
    @param f       a frequency grid
    @param t       a temperature profile

    @author Patrick Eriksson 08.04.2000 */
void planck (
              MATRIX&     B, 
        const VECTOR&     f,
        const VECTOR&     t );

/** Calculates the Planck function for a single temperature.

    @param B       output: the blackbody radiation
    @param f       a frequency grid
    @param t       a temperature value

    @author Patrick Eriksson 08.04.2000 */
void planck (
              VECTOR&     B, 
        const VECTOR&     f,
        const Numeric&    t );



//==========================================================================
//=== Tangent altitudes.
//==========================================================================

/** Calculates the geometrical tangent altitude.
    That is, refraction is neglected.

    @return        the tangent altitude
    @param view    the angle between zenith and the LOS
    @param z_plat  the platform altitude

    @author Patrick Eriksson 08.04.2000 */
Numeric ztan_geom(
        const Numeric&     view,
        const Numeric&     z_plat );



//==========================================================================
//=== Core functions for RTE and BL 
//==========================================================================

/** Performs a single iteration for RTE calculations (one viewing angle).
    The vector Y is not initilised, the obtained values are added to Y.
    Note that only a single iteration is performed. The ground is not 
    taken into account.     
    This function can be used to calculate emission spectra for parts of
    the atmosphere.
        
    @param y             output: the spectrum
    @param start_index   start index for the integration
    @param stop_index    stop index for the integration
    @param Tr            transmission matrix
    @param S             source function matrix
    @param n_f           number of frequencies

    @author Patrick Eriksson 15.06.2000 */
void rte_iterate (
             VECTOR&   y,
       const int&      start_index,
       const int&      stop_index,
       const MATRIX&   Tr,
       const MATRIX&   S,
       const size_t    n_f );



/** Performs the RTE calculations for one viewing angle.
    This function allows calculation of emission spectra for single
    viewing angles in function beside yRteXx.
        
    @param y             output: the spectrum
    @param start_index   start index for the integration
    @param stop_index    stop index for the integration
    @param Tr            transmission matrix
    @param S             source function matrix
    @param y_space       intensity entering the atmosphre at start of LOS
    @param ground        flag/index for ground intersection
    @param e_ground      ground emissivity
    @param y_ground      ground blackbody radiation 

    @author Patrick Eriksson 22.05.2000 */
void rte (
             VECTOR&   y,
       const int&      start_index,
       const int&      stop_index,
       const MATRIX&   Tr,
       const MATRIX&   S,
       const VECTOR&   y_space,
       const int&      ground,
       const VECTOR&   e_ground,
       const VECTOR&   y_ground );


/** Performs a single iteration for BL calculations (one viewing angle).
    The vector Y is not initilised, Y is multiplied with the obtained values.
    Note that only a single iteration is performed. The ground is not 
    taken into account.     
    This function can be used to calculate transmissions for parts of
    the atmosphere.
        
    @param y             output: the spectrum
    @param start_index   start index for the integration
    @param stop_index    stop index for the integration
    @param Tr            transmission matrix
    @param S             source function matrix
    @param n_f           number of frequencies

    @author Patrick Eriksson 15.06.2000 */
void bl_iterate (
             VECTOR&   y,
       const int&      start_index,
       const int&      stop_index,
       const MATRIX&   Tr,
       const size_t    n_f );



/** Performs the BL (transmission) calculations for one viewing angle.
    This function allows calculation of transmission spectra for single
    viewing angles in functions beside yBlXx.
        
    @param y             output: the spectrum
    @param start_index   start index for the integration
    @param stop_index    stop index for the integration
    @param Tr            transmission matrix
    @param ground        flag/index for ground intersection
    @param e_ground      ground emissivity

    @author Patrick Eriksson 22.05.2000 */
void bl (
             VECTOR&   y,
       const int&      start_index,
       const int&      stop_index,
       const MATRIX&   Tr,
       const int&      ground,
       const VECTOR&   e_ground );



//==========================================================================
//=== Conversion and interpolation of pressure and altitude grids.
//==========================================================================

/** Converts an altitude vector to pressures.
    The log. of the pressures are interpolated linearly.
    In Matlab notation:

      p = exp(interp1(z0,log(p0),z,'linear'))

    @param p       output: the pressures at z
    @param z0      original altitude grid
    @param p0      original pressure grid
    @param z       new altitude grid

    @author Patrick Eriksson 10.04.2000 */
void z2p(
              VECTOR&     p,
        const VECTOR&     z0,
        const VECTOR&     p0,
        const VECTOR&     z );

/** Interpolates a vertical profile at a new set of pressures.
    A linear interpolation using log. pressure is applied.
    In Matlab notation, the following expression is used:

      p = interp1(log(p0),x,log(p),'linear')

    @param x       output: the interpolated values at p
    @param p0      original pressure grid
    @param x0      the profile to be interpolated
    @param p       new pressure grid

    @author Patrick Eriksson 12.04.2000 */
void interpp(
              VECTOR&     x,
        const VECTOR&     p0,
        const VECTOR&     x0,
        const VECTOR&     p );

/** Interpolates a matrix, such as an absorption matrix, at a new 
    set of pressures.
    A linear interpolation using log. pressure is applied.
    In Matlab notation, the following expression is used:

      A = interp1(log(p0),A0,log(p),'linear')

    @param A       output: the interpolated values at p
    @param p0      original pressure grid
    @param A0      the matrix to be interpolated
    @param p       new pressure grid

    @author Patrick Eriksson 13.06.2000 */
void interpp(
              MATRIX&  A,
        const VECTOR&  p0, 
        const MATRIX&  A0, 
        const VECTOR&  p );
