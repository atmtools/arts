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



////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
   \file   m_wfs.cc

   This file contains functions to calculate WFs for a number of variables.

   Calculation of WFs are described in the AUG sections "Atmospheric WFs"
   and "Measurement errors".

   \author Patrick Eriksson
   \date 2000-09-14 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "arts.h"
#include "md.h"
#include "vecmat.h"
#include "messages.h"          
#include "wsv.h"          
#include "math_funcs.h"          
#include "atm_funcs.h"          
#include "hmatrix.h"
#include "los.h"
extern const Numeric PLANCK_CONST;
extern const Numeric BOLTZMAN_CONST;



////////////////////////////////////////////////////////////////////////////
//   Function(s) to join two WF matrices and related variables
////////////////////////////////////////////////////////////////////////////

//// k_append ////////////////////////////////////////////////////////////////
/**
   Appends the K matrix to either Kx or Kb.

   \retval   kx           the Kx matrix
   \retval   kx_names     identity names for KX
   \retval   kx_index     identity indecies KX
   \retval   kx_aux       additional data KX
   \param    k            the K matrix
   \param    k_names      identity names for K
   \param    k_aux        additional data K

   \author Patrick Eriksson
   \date   2000-10-20
*/
void k_append (
		    MATRIX&          kx,
		    ARRAYofstring&   kx_names,
		    ARRAYofsizet&    kx_lengths,
		    MATRIX&          kx_aux,
              const MATRIX&          k,
              const ARRAYofstring&   k_names,
              const MATRIX&          k_aux )
{
  // Size of Kx and K
  const size_t  ny1  = kx.nrows();         // length of measurement vector (y)
  const size_t  nx1  = kx.ncols();         // length of state vector (x)
  const size_t  nri1 = kx_names.size();    // number of retrieval identities
  const size_t  ny2  = k.nrows();  
  const size_t  nx2  = k.ncols();  
  const size_t  nri2 = k_names.size();
        size_t  iri;

         MATRIX   ktemp(ny1,nx1), ktemp_aux(nx1,2);
  ARRAYofstring   ktemp_names(nri1);
   ARRAYofsizet   ktemp_lengths(nri1);
  
  if ( nx1 > 0 )
  {
    if ( ny1 != ny2 )
      throw runtime_error(
            "The two WF matrices have different number of rows." ); 

    // Make copy of Kx data
    copy( kx,         ktemp );
    copy( kx_lengths, ktemp_lengths );
    copy( kx_aux,     ktemp_aux );
    copy( kx_names,   ktemp_names );
  }

  // Reallocate the Kx data
  resize( kx,         ny2,       nx1+nx2 );
  resize( kx_names,   nri1+nri2          );
  resize( kx_lengths, nri1+nri2          );
  resize( kx_aux,     nx1+nx2,   2       );

  // Move Ktemp to Kx
  if ( nx1 > 0 )
  {
    copy( ktemp,       kx.sub_matrix( 0, ny2, 0, nx1 ) );
    copy( ktemp_aux,   kx_aux.sub_matrix( 0, nx1, 0, 2 ) );
    for ( iri=0; iri<nri1; iri++ )
    {
      kx_lengths[iri]  = ktemp_lengths[iri];
      kx_names[iri]    = ktemp_names[iri];
    }
  }

  // Calculate the vector length for each identity in K
  size_t l = (size_t) floor(nx2/nri2);

  // Move K to Kx
  copy( k,       kx.sub_matrix( 0, ny2, nx1, nx1+nx2 ) );
  copy( k_aux,   kx_aux.sub_matrix( nx1, nx1+nx2, 0, 2 ) );
  for ( iri=0; iri<nri2; iri++ )
  {
    kx_names[nri1+iri]   = k_names[iri];
    kx_lengths[nri1+iri] = l;
  } 
}



////////////////////////////////////////////////////////////////////////////
//   Help functions for grid conversions
////////////////////////////////////////////////////////////////////////////

//// p2grid ////////////////////////////////////////////////////////////////
/**
   Converts a pressure grid to a grid fitting the grid functions below.

   -log(p) is used as altitude variable. The minus is included to get
    increasing values, a demand for the grid functions. 

   \retval   grid        -log(pgrid)
   \param    pgrid       pressure grid for some variable

   \author Patrick Eriksson
   \date   2001-01-15
*/
void p2grid(
	      VECTOR&   grid,
        const VECTOR&   pgrid )
{
  resize( grid, pgrid.size() );
  copy( pgrid, grid );
  transf( grid, log, grid );
  copy( scaled(grid,-1), grid );
}



//// grid2grid_index ///////////////////////////////////////////////////////
/**
   Gives indecies for conversion between vertical grids.

   This function gives the indeces of the grid points affected by a change at 
   a set of retrieval points. That is, a first step to calculate dkappa/dkp
   (see Equation 4 of the AUG section on "Basic atmospheric WFs). The second
   step is performed by the function GRID2GRID_WEIGHTS found below.
  
   This function can be used to convert LOS WFs to WFs for vertical altitudes.
   
   The profile/function for the retrieved quantity is treated to be piecewise
   linear. With other words, the base functions are tenth functions. The 
   functions areis treated to be 0 outside the end points. 

   The function returns a 2 column matrix, where column 1 is the lowest
   index affected and column 2 the highest index. If no points are affected
   -1 are returned. 
  
   An example:

   x0 = [1 2 3 4 5 6 7] and 

   xp = [0 1 2 3 4.5 5 5.5 6 8] gives

   Is = [-1    -1
          0     0
          1     1
          2     3
          3     3
          4     4
         -1    -1
          5     6
          6     6 ]

   Note that the XP point 5.5 gives -1 as this point not embraces any point 
   of X0 as there are retrieval points at the neighbouring X0 points.
   
   X0 is the original grid. For example, the LOS grid. 
   XP is the retrieval grid. For example, vertical altitudes.
  
   Both vectors most be sorted with increasing values. 
   The length of XP must be > 1.
  
   The algorithm used here is not the most effective. However, more elaborated
   algorithms were tested but they failed for some special case, and  instead 
   a simplistic approach was used.

   \retval   is          two column matrix with indecies of x0
   \param    x0          grid for some variable
   \param    xp          retrieval grid

   \author Patrick Eriksson
   \date   2000-09-15
*/
void grid2grid_index (
                    MATRIX&   is,
              const VECTOR&   x0,
              const VECTOR&   xp )
{
  const size_t n0 = x0.size();        // length of original grid
  const size_t np = xp.size();        // length if retrieval grid
        size_t i0, ip;                // counter for each grid

  // Resize is and set all values to -1
  resize( is, np, 2 );
  setto( is, -1.0 );
 
  // Some checks
  if ( np < 2 )
    throw logic_error("The retrieval grid must have a length > 1.");
  for ( i0=1; i0<n0; i0++ )
  {
    if ( x0[i0-1] >= x0[i0] )
      throw logic_error("The grids must be sorted with increasing values.");
  }
  for ( ip=1; ip<np; ip++ )
  {
    if ( xp[ip-1] >= xp[ip] )
      throw logic_error("The grids must be sorted with increasing values.");
  }

  // Do something only if at least one point of X0 is inside XP
  if ( !( (x0[0]>xp[np-1]) || (x0[n0-1]<xp[0]) ) )
  {

    // First point of XP
    i0 = 0;
    for ( ; x0[i0] < xp[0]; i0++ ) {}
    if ( x0[i0] < xp[1] )
    {
      is[0][0] = (Numeric) i0;
      for ( ; (i0<n0) && (x0[i0]<xp[1]); i0++ ) {}
      is[0][1] = (Numeric) (i0 - 1);
    }

    // Points inside XP
    for ( ip=1; ip<(np-1); ip++ )
    {
      i0 = 0;
      for ( ; (i0<n0) && (x0[i0]<=xp[ip-1]); i0++ ) {}
      if ( (i0<n0) && (x0[i0]<xp[ip+1]) )
      {
        is[ip][0] = (Numeric) i0;
        for ( ; (i0<n0) && (x0[i0]<xp[ip+1]); i0++ ) {}
        is[ip][1] = (Numeric) (i0 - 1);
      }
    }

    // Last point of XP
    i0 = 0;
    for ( ; (i0<n0) && (x0[i0]<=xp[np-2]); i0++ ) {}
    if ( (i0<n0) && (x0[i0]<=xp[np-1]) )
    {
      is[np-1][0] = (Numeric) i0;
      for ( ; (i0<n0) && (x0[i0]<=xp[np-1]); i0++ ) {}
      is[np-1][1] = (Numeric) (i0 - 1);
    }
  }
}



//// grid2grid_weights /////////////////////////////////////////////////////
/**
   Gives weights for conversion between vertical grids.
  
   This function returns the derivative dkappa/dkp (see the AUG section on 
   "Atmospheric WFs"), here dennoted as the weight, for the indeces selected.

   This is an accompanying function to GRID2GRID_INDEX. This latter function
   gives the indeces for which the weight > 0.
  
   The weight is

   w = abs( (x0-xp(ip+-1)) / (xp(ip)-xp(ip+-1)))

   where x0 is the point of the original grid (kappa), xp(ip) is the 
   considered retrieval point and xp(ip+-1) is the retrieval point on the
   other side of x0.
  
   \retval   w    weights for each point in x0
   \param    x0   the original grid, e.g. the LOS grid. 
   \param    i1   first index of x0 to consider
   \param    i2   last index of x0 to consider
   \param    xp   retrieval grid, e.g. vertical altitudes
   \param    ip   index of retrieval point

   \author Patrick Eriksson
   \date   2000-09-15
*/
void grid2grid_weights (
                    VECTOR&   w,
              const VECTOR&   x0,
              const size_t&   i1,
              const size_t&   i2,
              const VECTOR&   xp,
              const size_t&   ip )
{
  const size_t   np = xp.size();        // number of retrieval points
  const size_t   nw = i2-i1+1;          // number of LOS points affected
        size_t   i;

  // Reallocate w
  resize(w,nw);

  // First point of the retrieval grid
  if ( ip == 0 )
  {
    for ( i=0; i<nw; i++ )
      w[i] = ( xp[1] - x0[i1+i] ) / ( xp[1] - xp[0] );
  }
  // Points inside the retrieval grid
  else if ( ip < (np-1) )  
  {
    for ( i=0; i<nw; i++ )
    {
      if ( x0[i1+i] <= xp[ip] )
        w[i] = ( x0[i1+i] - xp[ip-1] ) / ( xp[ip] - xp[ip-1] );
      else
        w[i] = ( xp[ip+1] - x0[i1+i] ) / ( xp[ip+1] - xp[ip] );
    }
  }
  // Last point of the retrieval grid
  else
  {
    for ( i=0; i<nw; i++ )
      w[i] = ( x0[i1+i] - xp[np-2] ) / ( xp [np-1] - xp[np-2] );
  }
}



////////////////////////////////////////////////////////////////////////////
//   Help functions for absloswfsCalc to calculate absorption LOS WFs
////////////////////////////////////////////////////////////////////////////

//// absloswfs_rte_1pass /////////////////////////////////////////////////////
/**
   Calculates absorption LOS WFs f single pass cases.

   Help function for absloswfsCalc treating a single zenith angle.

   This function covers cases where each point is passed once.
   That is 1D limb sounding and downward observatiuons are not covered.

   The expression used are described in sub-section 2.1 of the AUG section
   "Atmospheric WFs".

   \retval   k             abs. LOS WFs
   \param    y             spectrum vector
   \param    start_index   start LOS index for iteration
   \param    stop_index    stop LOS index for iteration
   \param    lstep         length between LOS points
   \param    tr            transmissions
   \param    s             source function values
   \param    ground        ground flag
   \param    e_ground      ground emissivity
   \param    y_ground      ground emission

   \author Patrick Eriksson
   \date   2000-09-15
*/
void absloswfs_rte_1pass (
                    MATRIX&   k,
                    VECTOR    y,
              const size_t&   start_index,
	      const size_t&   stop_index,   // this variable is used by
              const Numeric&  lstep,        // absloswfs_down
              const MATRIX&   tr,
              const MATRIX&   s,
              const INDEX&    ground,
              const VECTOR&   e_ground,
              const VECTOR&   y_ground )
{
  const size_t   nf = tr.nrows();     // number of frequencies
        size_t   iv;                  // frequency index
        VECTOR   t(nf,1.0);           // transmission to the sensor
        size_t   q;                   // LOS index (same notation as in AUG)

  if ( ground && ((ground-1==stop_index) || (ground-1==start_index)) )
    throw logic_error("The ground cannot be at one of the end points."); 

  // Resize K
  resize( k, nf, start_index+1 );

  // We start here at the LOS point closest to the sensor, that is,
  // reversed order compared to RTE_ITERATE  

  // The LOS point closest to the sensor 
  q  = stop_index;
  for ( iv=0; iv<nf; iv++ )    
  {
    t[iv]   *= tr[iv][q];
    y[iv]    = (y[iv]-s[iv][q]*(1.0-tr[iv][q]))/tr[iv][q];
    k[iv][q] = (-lstep/2)*(y[iv]-s[iv][q])*t[iv];
  }

  // Points inside the LOS
  for ( q=stop_index+1; q<start_index; q++ )
  {
    // Not a ground point
    if ( q != (ground-1) )
    {
      for ( iv=0; iv<nf; iv++ )    
      {
        y[iv]    = (y[iv]-s[iv][q]*(1.0-tr[iv][q]))/tr[iv][q];
        k[iv][q] = (-lstep/2)*( 2*(y[iv]-s[iv][q])*tr[iv][q] + 
                                             s[iv][q] - s[iv][q-1] ) *t[iv];
        t[iv]   *= tr[iv][q]; 
      }
    }
    // A ground point
    else
    {
      out1 << "WARNING: The function absloswfs_1pass not tested for "
              "ground reflections\n";
      for ( iv=0; iv<nf; iv++ )    
      {
        y[iv]    = ( y[iv] - e_ground[iv]*y_ground[iv] - 
                            s[iv][q]*(1.0-tr[iv][q])*(1-e_ground[iv]) ) / 
                                             ( tr[iv][q] * (1-e_ground[iv]) );
        k[iv][q] = (-lstep/2)*( 2*(y[iv]-s[iv][q])*tr[iv][q]*(1-e_ground[iv])+ 
                                s[iv][q]*(1-e_ground[iv]) + 
                             e_ground[iv]*y_ground[iv] - s[iv][q-1] ) * t[iv];
        t[iv]   *= tr[iv][q] * (1-e_ground[iv]); 
      }
    }
  }

  // The LOS point furthest away from the sensor
  q = start_index;
  for ( iv=0; iv<nf; iv++ )    
    k[iv][q]  = (-lstep/2)*(y[iv]-s[iv][q-1])*t[iv];

  // To check the function: Y shall now be equal to the radiation entering 
  // the atmosphere and T the total transmission
}



//// absloswfs_rte_limb ///////////////////////////////////////////////////////
/**
   Calculates absorption LOS WFs for 1D limb sounding.

   Help function for absloswfsCalc treating a single zenith angle.

   The expression used are described in sub-section 2.2 of the AUG section
   "Atmospheric WFs".

   \retval   k             abs. LOS WFs
   \param    y             spectrum vector
   \param    y_space       cosmic radiation
   \param    start_index   start LOS index for iteration
   \param    lstep         length between LOS points
   \param    tr            transmissions
   \param    s             source function values
   \param    ground        ground flag
   \param    e_ground      ground emissivity

   \author Patrick Eriksson
   \date   2000-09-15
*/
void absloswfs_rte_limb (
                    MATRIX&   k,
                    VECTOR    y,              // = y_q^q
              const VECTOR&   y_space,             
              const size_t&   start_index,
              const Numeric&  lstep,
              const MATRIX&   tr,
              const MATRIX&   s,
              const INDEX&    ground,
              const VECTOR&   e_ground )
{
  const size_t   nf = tr.nrows();     // number of frequencies
        size_t   iv;                  // frequency index
        VECTOR   t1q;                 // transmission tangent point - q squared
        VECTOR   tqn(nf,1.0);         // transmission q - sensor
        size_t   q;                   // LOS index (same notation as in AUG)
       Numeric   tv, tv1;             // transmission value for q and q-1
        VECTOR   yn(nf);              // = y_space

  copy( y_space, yn );

  // Resize K
  resize( k, nf, start_index+1 );

  // Calculate the total transmission
  bl( t1q, start_index, start_index, tr, ground, e_ground );

  // We start at the outermost point
  q  = start_index;       
  for ( iv=0; iv<nf; iv++ )    
  {
    tv1      = tr[iv][q-1];
    t1q[iv] /= tv1*tv1;
    tqn[iv] *= tv1;
    y[iv]    = ( y[iv] - s[iv][q-1]*(1-tv1)*(1+t1q[iv]*tv1) ) / tv1;
    k[iv][q]  = (-lstep/2)*( ( 2*yn[iv] + s[iv][q-1]*(1-2*tv1) ) * 
                              t1q[iv]*tv1 + y[iv] - s[iv][q-1] )*tv1;
  }

  // Points inside the LOS
  for ( q=start_index-1; q>0; q-- )
  {
    for ( iv=0; iv<nf; iv++ )    
    {
      tv1      = tr[iv][q-1];    
      tv       = tr[iv][q];
      t1q[iv] /= tv1*tv1;
      y[iv]    = ( y[iv] - s[iv][q-1]*(1-tv1)*(1+t1q[iv]*tv1) ) / tv1;
      k[iv][q] = (-lstep/2) * ( 
           ( 4*(yn[iv]-s[iv][q])*tv1*tv + 3*(s[iv][q]-s[iv][q-1])*tv1 + 
                                        2*s[iv][q-1] ) * t1q[iv]*tv1 + 
             2*(y[iv]-s[iv][q-1])*tv1 + s[iv][q-1] - s[iv][q] ) * tqn[iv];
      tqn[iv] *= tv1;
      yn[iv]   = yn[iv]*tv + s[iv][q]*(1-tv);
    } 
  }

  // The tangent or ground point
  for ( iv=0; iv<nf; iv++ )    
    k[iv][0]  = (-lstep/2)*( (2*yn[iv]*tv1+s[iv][0]*(1-2*tv1))*t1q[iv] +
                       y[iv] - s[iv][0] ) * tqn[iv];

  // To check the function
  // Without ground reflection: T1Q=1 and Y=0
  // With ground reflection: T1Q=(1-e) and Y=eB
  // print_vector( t1q );
  // print_vector( y );
}



//// absloswfs_rte_down //////////////////////////////////////////////////////
/**
   Calculates absorption LOS WFs for 1D downward looking observations.

   Help function for absloswfsCalc treating a single zenith angle.

   The expression used are described in sub-section 2.3 of the AUG section
   "Atmospheric WFs".

   \retval   k             abs. LOS WFs
   \param    y             spectrum vector
   \param    y_space       cosmic radiation
   \param    start_index   start LOS index for iteration
   \param    stop_index    stop LOS index for iteration
   \param    lstep         length between LOS points
   \param    tr            transmissions
   \param    s             source function values
   \param    ground        ground flag
   \param    e_ground      ground emissivity
   \param    y_ground      ground emission

   \author Patrick Eriksson
   \date   2000-09-15
*/
void absloswfs_rte_down (
                    MATRIX&   k,
              const VECTOR&   y,
              const VECTOR&   y_space,
              const size_t&   start_index,
              const size_t&   stop_index,
              const Numeric&  lstep,
              const MATRIX&   tr,
              const MATRIX&   s,
              const INDEX&    ground,
              const VECTOR&   e_ground,
              const VECTOR&   y_ground )
{
  const size_t   nf = tr.nrows(); // number of frequencies
        size_t   iv;              // frequency index
        size_t   q;               // LOS index (same notation as in AUG)
        VECTOR   y0(nf);          // see below
        MATRIX   k2;              // matrix for calling other LOS WFs functions
        VECTOR   tr0;             // see below

  // Resize K
  resize( k, nf, start_index+1 );

  // Calculate Y0, the intensity reaching the platform altitude at the far
  // end of LOS, that is, from above.
  copy( y_space, y0 );
  rte_iterate( y0, start_index-1, stop_index, tr, s, nf );

  // Calculate TR0, the transmission from the platform altitude down to the
  // tangent point or the ground, and up to the platform again.
  bl( tr0, stop_index, stop_index, tr, ground, e_ground );

  // The indeces below STOP_INDEX are handled by the limb sounding function.
  // The limb function is given Y0 instead of cosmic radiation 
  absloswfs_rte_limb( k2, y, y0, stop_index, lstep, tr, s, ground, e_ground );
  for ( iv=0; iv<nf; iv++ )
  {
    for ( q=0; q<stop_index; q++ )
      k[iv][q] = k2[iv][q];
  }

  // The indeces above STOP_INDEX are handled by the 1pass function.
  // The transmission below STOP_INDEX must here be considered.
  absloswfs_rte_1pass( k2, y0, start_index, stop_index, lstep, tr, s, 
                     ground, e_ground, y_ground );
  for ( iv=0; iv<nf; iv++ )
  {
    for ( q=stop_index+1; q<=start_index; q++ )
      k[iv][q] = k2[iv][q]*tr0[iv];
  }

  // The platform altitude must be treated seperately
  //
  // Calculate the intensity generated below point q-1, YQQ
  VECTOR yqq(nf);
  rte( yqq, stop_index-1, stop_index-1, tr, s, VECTOR(nf,0.0), 
                                            ground, e_ground, y_ground );
  //
  // Y0 is moved one step upwards and TR0 one step downwards
  q = stop_index; 
  for ( iv=0; iv<nf; iv++ )
  {
    tr0[iv] /= tr[iv][q-1]*tr[iv][q-1];    
    y0[iv]   = ( y0[iv] - s[iv][q]*(1-tr[iv][q]) ) / tr[iv][q];
    k[iv][q] = (-lstep/2)*( (3*(y0[iv]-s[iv][q])*tr[iv][q-1]*tr[iv][q] + 
                  2*(s[iv][q]-s[iv][q-1])*tr[iv][q-1] + s[iv][q-1] )*tr0[iv] + 
                                         yqq[iv] - s[iv][q-1] )*tr[iv][q-1];
  }
}



//// absloswfs_rte ///////////////////////////////////////////////////////////
/**
   \author Patrick Eriksson
   \date   2000-09-15
*/
void absloswfs_rte (
                    ARRAYofMATRIX&   absloswfs,
              const LOS&             los,   
              const ARRAYofMATRIX&   source,
              const ARRAYofMATRIX&   trans,
              const VECTOR&          y,
              const VECTOR&          y_space,
              const VECTOR&          f_mono,
              const VECTOR&          e_ground,
              const Numeric&         t_ground )
{
  const size_t  nza = los.start.size();   // number of zenith angles  
  const size_t  nf  = f_mono.size();      // number of frequencies  
        VECTOR  yp(nf);                   // part of Y
        size_t  iy0=0;                    // y index

  // Set up vector for ground blackbody radiation
  VECTOR   y_ground(f_mono.size()); 
  if ( any_ground(los.ground) )
    planck( y_ground, f_mono, t_ground );

  // Resize the LOS WFs array
  resize( absloswfs, nza );

  // Loop zenith angles
  out3 << "    Zenith angle nr:\n      ";
  for (size_t i=0; i<nza; i++ ) 
  {
    if ( ((i+1)%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    // Do something only if LOS has any points
    if ( los.p[i].size() > 0 )
    {

      // Pick out the part of Y corresponding to the present zenith angle
      copy( y(iy0,iy0+nf), yp );

      // The calculations are performed in 3 sub-functions
      //
      // Upward looking (=single pass)
      if ( los.stop[i]==0 )
        absloswfs_rte_1pass( absloswfs[i], yp, los.start[i], 0, los.l_step[i], 
                     trans[i], source[i], los.ground[i], e_ground, y_ground );

      //
      // 1D limb sounding
      else if ( los.start[i] == los.stop[i] )
        absloswfs_rte_limb( absloswfs[i], yp, y_space, los.start[i], 
                los.l_step[i], trans[i], source[i], los.ground[i], e_ground );

      //
      // 1D downward looking
      else 
        absloswfs_rte_down( absloswfs[i], yp, y_space, los.start[i], 
                            los.stop[i], los.l_step[i], trans[i], source[i], 
                                           los.ground[i], e_ground, y_ground );
    }

    iy0 += nf;        
  }

  out3 << "\n";
}



//// absloswfs_tau ///////////////////////////////////////////////////////////
/**
   \author Patrick Eriksson
   \date   2001-03-30
*/
void absloswfs_tau (
                    ARRAYofMATRIX&   absloswfs,
	      const LOS&             los,
	      const VECTOR&          f_mono )
{
  const size_t  nza = los.start.size();   // number of zenith angles  
  const size_t  nf  = f_mono.size();      // number of frequencies  
        Numeric kw, kw2;
        INDEX   row, col, np;

  // Resize the LOS WFs array
  resize( absloswfs, nza );

  // Loop zenith angles
  out3 << "    Zenith angle nr:\n      ";
  for (size_t i=0; i<nza; i++ ) 
  {
    if ( ((i+1)%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();

    np = los.start[i] + 1;

    resize( absloswfs[i], nf, np );
    
    // Do something only if LOS has any points
    if ( los.p[i].size() > 0 )
    {

      // Upward looking (=single pass)
      if ( los.stop[i]==0 )
      {
        kw  = los.l_step[i] / 2.0;
        kw2 = los.l_step[i];
        for( row=0; row<nf; row++ )
	{
          absloswfs[i][row][0] = kw;
          absloswfs[i][row][np-1] = kw;
          for( col=1; col<np-1; col++ )
            absloswfs[i][row][col] = kw2;
        }
      }

      // 1D limb sounding
      else if ( los.start[i] == los.stop[i] )
      {
        kw  = los.l_step[i];
        kw2 = los.l_step[i] * 2.0;
        for( row=0; row<nf; row++ )
	{
          absloswfs[i][row][0] = kw;
          absloswfs[i][row][np-1] = kw;
          for( col=1; col<np-1; col++ )
            absloswfs[i][row][col] = kw2;
        }
      }

      // 1D downward looking
      else 
      {
        kw  = los.l_step[i];
        kw2 = los.l_step[i] * 2.0;
        for( row=0; row<nf; row++ )
	{
          absloswfs[i][row][0] = kw;
          absloswfs[i][row][np-1] = kw / 2.0;
          absloswfs[i][row][los.stop[i]] = kw * 1.5;
          for( col=1; col<los.stop[i]; col++ )
            absloswfs[i][row][col] = kw2;
          for( col=los.stop[i]+1; col<np-1; col++ )
            absloswfs[i][row][col] = kw;
        }
      }
    }
  }
  out3 << "\n";
}



////////////////////////////////////////////////////////////////////////////
//   Functions to calculate source function LOS WFs
////////////////////////////////////////////////////////////////////////////

//// sourceloswfs_1pass ////////////////////////////////////////////////////
/**
   Calculates source function LOS WFs for single pass cases.

   Help function for sourceloswfs treating a single zenith angle.

   This function covers cases where each point is passed once.
   That is 1D limb sounding and downward observatiuons are not covered.

   The expression used are described in sub-section 3.1 of the AUG section
   "Atmospheric WFs".

   \retval   k             source LOS WFs
   \param    start_index   start LOS index for iteration
   \param    stop_index    stop LOS index for iteration
   \param    tr            transmissions
   \param    ground        ground flag
   \param    e_ground      ground emissivity

   \author Patrick Eriksson
   \date   2000-09-15
*/
void sourceloswfs_1pass (
                    MATRIX&   k,
              const size_t&   start_index,
	      const size_t&   stop_index,   // this variable is used by 1D down
              const MATRIX&   tr,
              const INDEX&    ground,
              const VECTOR&   e_ground )
{
  const size_t   nf = tr.nrows();     // number of frequencies
        size_t   iv;                  // frequency index
        VECTOR   t(nf,1.0);           // transmission to the sensor
        size_t   q;                   // LOS index (same notation as in AUG) 

  if ( ground && ((ground-1==stop_index) || (ground-1==start_index)) )
    throw logic_error("The ground cannot be at one of the end points."); 

  // Resize K
  resize( k, nf, start_index+1 );

  // We start here at the LOS point closest to the sensor, that is,
  // reversed order compared to RTE_ITERATE  

  // The LOS point closest to the sensor 
  q  = stop_index;
  for ( iv=0; iv<nf; iv++ )    
    k[iv][q]  = ( 1 - tr[iv][q] ) / 2;

  // Points inside the LOS
  for ( q=stop_index+1; q<start_index; q++ )
  {
    // Not a ground point
    if ( q != ground )
    {
      for ( iv=0; iv<nf; iv++ )    
      {
        k[iv][q] = ( 1 - tr[iv][q-1]*tr[iv][q] ) * t[iv] / 2;
        t[iv]  *= tr[iv][q]; 
      }
    }
    // A ground point
    else
    {
      out1 << "WARNING: The function sourceloswfs_1pass not tested "
              "for ground reflections\n";

      for ( iv=0; iv<nf; iv++ )    
      {
        k[iv][q] = ( (1-tr[iv][q])*(1-e_ground[iv])*tr[iv][q-1] + 1 - 
                                                     tr[iv][q-1] ) * t[iv] / 2;
        t[iv]  *= tr[iv][q]*(1-e_ground[iv]); 
      }
    }
  }

  // The LOS point furthest away from the sensor
  q = start_index;
  for ( iv=0; iv<nf; iv++ )    
    k[iv][q]  = ( 1 - tr[iv][q-1] ) * t[iv] / 2;
}



//// sourceloswfs_limb /////////////////////////////////////////////////////
/**
   Calculates source function LOS WFs for 1D limb sounding.

   Help function for sourceloswfs treating a single zenith angle.

   The expression used are described in sub-section 3.2 of the AUG section
   "Atmospheric WFs".

   \retval   k             abs. LOS WFs
   \param    start_index   start LOS index for iteration
   \param    tr            transmissions
   \param    ground        ground flag
   \param    e_ground      ground emissivity

   \author Patrick Eriksson
   \date   2000-09-15
*/
void sourceloswfs_limb (
                    MATRIX&   k,
              const size_t&   start_index,
              const MATRIX&   tr,
              const INDEX&    ground,
              const VECTOR&   e_ground )
{
  const size_t   nf = tr.nrows();      // number of frequencies
        size_t   iv;                  // frequency index
        VECTOR   t1q;                 // transmission tangent point - q squared
        VECTOR   tqn(nf,1);           // transmission q - sensor
        size_t   q;                   // LOS index (same notation as in AUG)

  // Calculate the total transmission
  bl( t1q, start_index, start_index, tr, ground, e_ground );

  // Resize K
  resize( k, nf, start_index+1 );

  // We start at the outermost point
  q  = start_index;       
  for ( iv=0; iv<nf; iv++ )    
  {
    t1q[iv] /= tr[iv][q-1] * tr[iv][q-1];
    k[iv][q]  = ( (1-tr[iv][q-1])*t1q[iv]*tr[iv][q-1] + 1 - tr[iv][q-1] )/2;
  }

  // Points inside the LOS
  for ( q=start_index-1; q>0; q-- )
  {
    for ( iv=0; iv<nf; iv++ )    
    {
      t1q[iv]  /= tr[iv][q-1] * tr[iv][q-1];
      k[iv][q]  = ( (1-tr[iv][q-1]*tr[iv][q])*t1q[iv]*tr[iv][q-1]*
         tr[iv][q] + (1-tr[iv][q-1])*tr[iv][q] + 1 - tr[iv][q] ) * tqn[iv] / 2;
      tqn[iv]  *= tr[iv][q-1];
    } 
  }

  // The tangent or ground point
  for ( iv=0; iv<nf; iv++ )    
    k[iv][0]  = ( (1-tr[iv][0])*(1+t1q[iv]*tr[iv][0]) ) * tqn[iv] / 2;
}



//// sourceloswfs_down //////////////////////////////////////////////////////
/**
   Calculates source function LOS WFs for 1D downward looking observations.

   Help function for sourceloswfs treating a single zenith angle.

   The expression used are described in sub-section 3.3 of the AUG section
   "Atmospheric WFs".

   \retval   k             abs. LOS WFs
   \param    start_index   start LOS index for iteration
   \param    stop_index    stop LOS index for iteration
   \param    tr            transmissions
   \param    ground        ground flag
   \param    e_ground      ground emissivity

   \author Patrick Eriksson
   \date   2000-09-15
*/
void sourceloswfs_down (
                    MATRIX&   k,
              const size_t&   start_index,
              const size_t&   stop_index,
              const MATRIX&   tr,
              const INDEX&    ground,
              const VECTOR&   e_ground )
{
  const size_t   nf = tr.nrows(); // number of frequencies
        size_t   iv;             // frequency index
        size_t   q;              // LOS index (same notation as in AUG)
        MATRIX   k2;             // matrix for calling other LOS WFs functions
        VECTOR   tr0;            // see below

  // Resize K
  resize( k, nf, start_index+1 );

  // Calculate TR0, the transmission from the platform altitude down to the
  // tangent point or the ground, and up to the platform again.
  bl( tr0, stop_index, stop_index, tr, ground, e_ground );

  // The indeces below STOP_INDEX are handled by the limb sounding function.
  sourceloswfs_limb( k2, stop_index, tr, ground, e_ground );  
  for ( iv=0; iv<nf; iv++ )
  {
    for ( q=0; q<stop_index; q++ )
      k[iv][q] = k2[iv][q];
  }

  // The indecies above STOP_INDEX are handled by the 1pass function.
  // The transmission below STOP_INDEX must here be considered.
  sourceloswfs_1pass( k2, start_index, stop_index, tr, ground, e_ground );
  for ( iv=0; iv<nf; iv++ )
  {
    for ( q=stop_index+1; q<=start_index; q++ )
      k[iv][q] = k2[iv][q]*tr0[iv];
  }

  // The platform altitude must be treated seperately
  //
  // TR0 is moved one step downwards
  q = stop_index; 
  for ( iv=0; iv<nf; iv++ )
  {
    tr0[iv] /= tr[iv][q-1]*tr[iv][q-1];    
    k[iv][q]  = ( (1-tr[iv][q-1]*tr[iv][q])*tr0[iv]*tr0[iv]*tr[iv][q-1] + 1 - 
                 tr[iv][q-1] ) / 2;
  }
}



//// sourceloswfs //////////////////////////////////////////////////////////
/**
   Calculates source function LOS WFs,

   Main function for calulation of source function LOS WFs..

   The expression used are described in sub-section 2.3 of the AUG section
   "Atmospheric WFs".

   \retval   sourceloswfs   source LOS WFs
   \param    los            line of sight structure
   \param    tr             transmissions
   \param    e_ground       ground emissivity

   \author Patrick Eriksson
   \date   2000-09-15
*/
void sourceloswfs (
                    ARRAYofMATRIX&   sourceloswfs,
              const LOS&             los,   
              const ARRAYofMATRIX&   trans,
              const VECTOR&          f_mono,
              const VECTOR&          e_ground )
{
  const size_t  nza = los.start.size();   // number of zenith angles  

  // Resize the LOS WFs array
  resize(sourceloswfs,nza);

  // Loop zenith angles
  out3 << "    Zenith angle nr:      ";
  for (size_t i=0; i<nza; i++ ) 
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    // Do something only if LOS has any points
    if ( los.p[i].size() > 0 )
    {
      // The calculations are performed in 3 sub-functions
      //
      // Upward looking (=single pass)
      if ( los.stop[i]==0 )
        sourceloswfs_1pass( sourceloswfs[i], los.start[i], 0, trans[i], 
                                                     los.ground[i], e_ground );
      //
      // 1D limb sounding
      else if ( los.start[i] == los.stop[i] )
        sourceloswfs_limb( sourceloswfs[i], los.start[i], trans[i], 
                                                     los.ground[i], e_ground );
      //
      // 1D downward looking
      else 
        sourceloswfs_down( sourceloswfs[i], los.start[i], los.stop[i], 
                                           trans[i], los.ground[i], e_ground );
    }
  }
  out3 << "\n";
}



////////////////////////////////////////////////////////////////////////////
//   Core functions for analytical WFs.
//     Analytical expressions are used for species, temperature without 
//     hydrostatic eq. and continuum absorption.
//     These functions are very similar and a change in one of the functions
//     should most likely be included in the other functions.
////////////////////////////////////////////////////////////////////////////

//// k_species /////////////////////////////////////////////////////////////
/**
   Species WFs for one or several species.

   Calculates semi-analytical WFs for one or several species using 
   precalculated absorption LOS WFs. 

   The expression used are described in sub-section 5 of the AUG section
   "Atmospheric WFs".

   The avaliable units are
     1 fractions of linearisation state 
     2 volume mixing ratio
     3 number density

   \retval   k            weighting function matrix
   \retval   k_names      identity name(s)
   \retval   k_aux        additional data
   \param    los          line of sight structure
   \param    absloswfs    absorption LOS Wfs
   \param    p_abs        pressure grid for abs. calculations
   \param    t_abs        temperatures at p_abs
   \param    tags         all tags
   \param    abs_per_tg   absorption for each tag
   \param    vmrs         VMR profiles at p_abs
   \param    k_grid       retrieval grid
   \param    tg_nr        WFs are calculated for these tag index
   \param    unit         unit for the WFs (see above)   

   \author Patrick Eriksson
   \date   2000-09-15

    Adapted to MTL. 
    \date   2001-01-05
    \author Stefan Buehler
*/
void k_species (
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_aux,
              const LOS&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          p_abs,
              const VECTOR&          t_abs,             
              const TagGroups&       tags,
              const ARRAYofMATRIX&   abs_per_tg,
              const ARRAYofVECTOR&   vmrs,
              const VECTOR&          k_grid,
              const ARRAYofsizet&    tg_nr,
              const string&          unit )
{
  // Main sizes
  const size_t  nza = los.start.size();      // number of zenith angles  
  const size_t  nv  = abs_per_tg[0].nrows(); // number of frequencies
  const size_t  ntg = tg_nr.size();          // number of retrieval tags to do
  const size_t  np  = k_grid.size();         // number of retrieval altitudes

  // -log(p) is used as altitude variable. The minus is included to get
  // increasing values, a demand for the grid functions. 
  VECTOR  lgrid;                      // -log of the retrieval pressures
  p2grid( lgrid, k_grid );
  VECTOR  lplos;                      // -log of the LOS pressures

  // Indices
  // IP0 and IF0 are the index off-sets for the total K matrix
        size_t  itg;                       // Tag index
        size_t  iza;                       // Zenith angle index
        size_t  ip, ip0=0;                 // Retrieval point indices
        size_t  iv, iv0;                   // Frequency indices
        size_t  i1, iw;                    // weight indices

  // Other variables
        MATRIX  abs;                       // absorption at k_grid
        MATRIX  is;                        // matrix for storing LOS index
        VECTOR  w;                         // weights for LOS WFs
        VECTOR  a(nv);                     // temporary vector
        size_t  tg;                        // present tag nr
        VECTOR  vmr, p, t;                 // for conversion to VMR and ND 
        VECTOR  nd;                        // number density

  // Set up K and additional data. Set all values of K to 0
  resize(k,nza*nv,ntg*np);
  setto( k, 0.0);
  resize(k_names,ntg);
  resize(k_aux,ntg*np,2);

  // The calculations
  // Loop order:
  //   1 tags
  //   2 zenith angle
  //   3 retrieval altitudes 
  //   4 frequencies
  for ( itg=0; itg<ntg; itg++ ) 
  {
    // Present tag nr
    tg = tg_nr[itg];

    out2 << "  Doing " << tags[tg][0].Name() << "\n";

    // Get the absorption of the species at the retrieval points
    resize( abs, nv, np );
    interpp( abs, p_abs, abs_per_tg[tg], k_grid );

    // Fill K_NAMES
    {
      ostringstream os;
      os << "Species: " << tags[tg][0].Name();
      k_names[itg] = os.str();
    }

    // Fill column 0 of K_AUX
    for ( ip=0; ip<np; ip++ )
       k_aux[ip0+ip][0] = k_grid[ip];

    // Fill column 1 of K_AUX
    //   frac : fractions of the linearisation profile
    //   vmr  : VMR
    //   nd   : number density
    //
    resize( vmr, np );
    interpp( vmr, p_abs, vmrs[tg], k_grid ); 
    //
    if ( unit == "frac" )
    {
      for ( ip=0; ip<np; ip++ )
        k_aux[ip0+ip][1] = 1.0;
    }
    else if ( unit == "vmr" )
    {
      for ( ip=0; ip<np; ip++ )
      {
        for ( iv=0; iv<nv; iv++ )
          abs[iv][ip]   /= vmr[ip];
        k_aux[ip0+ip][1] = vmr[ip];
      }
    }  
    else if ( unit == "nd" )
    {
      resize( nd, np );
      interpp(  nd, p_abs, number_density(p_abs,t_abs), k_grid );
      ele_mult( vmr, nd, nd );
      for ( ip=0; ip<np; ip++ )
      {
        for ( iv=0; iv<nv; iv++ )
          abs[iv][ip]   /= nd[ip];
        k_aux[ip0+ip][1] = nd[ip];
      }
    }
    else
      throw runtime_error(
        "Allowed retrieval units are \"frac\", \"vmr\" and \"nd\"."); 

    // Set frequency zenith angle index off-set to 0
    iv0 = 0;                 

    // Loop zenith angles
    out3 << "    Zenith angle nr:\n      ";
    for ( iza=0; iza<nza; iza++ ) 
    {
      if ( ((iza+1)%20)==0 )
        out3 << "\n      ";
      out3 << " " << iza; cout.flush();
      
      // Do something only if LOS has any points
      if ( los.p[iza].size() > 0 )
      {
        // Get the LOS points affected by each retrieval point
	//        lplos = -1.0 * log(los.p[iza]);
        p2grid( lplos, los.p[iza] );
        grid2grid_index( is, lplos, lgrid );

        // Loop retrieval points
        for ( ip=0; ip<np; ip++ ) 
        {
          // Check if there is anything to do
          if ( is[ip][0] >= 0 )
          {
            // Get the weights for the LOS points
	    grid2grid_weights( w, lplos, (size_t)is[ip][0], (size_t)is[ip][1], 
                                                                 lgrid, ip );

            // Calculate the WFs.
            // A is di/dkappa*dkappa/dkp in a compact form.
            // This is possible as the columns of dkappa/dkp are identical.  
            setto( a, 0.0 );                     
            i1 = (size_t)is[ip][0];       // first LOS point to consider
            for ( iv=0; iv<nv; iv++ )
	    {
              for ( iw=i1; iw<=(size_t)is[ip][1]; iw++ )
                a[iv] += absloswfs[iza][iv][iw] * w[iw-i1];
              k[iv0+iv][ip0+ip] = a[iv] * abs[iv][ip];                    
	    }
          }            
        }
      }

      // Update the frequency index offset
      iv0 += nv;
    }  
    out3 << "\n";

    // Increase retrieval altitude index off-set 
    ip0 += np;
  }  
}



//// k_contabs //////////////////////////////////////////////////////////////
/**
   WFs for fit of continuum absorption.

   Calculates semi-analytical WFs for fit of continuum absorption using
   precalculated absorption LOS WFs. 

   The expression used are described in sub-section 6 of the AUG section
   "Atmospheric WFs".

   The continuum is fitted be determining an off-set at a number of
   points (order+1) that are evenly spread between the frequencies
   specified by the two index given. The weighting functions are zero
   outside these limits.

   \retval   k            weighting function matrix
   \retval   k_names      identity name(s)
   \retval   k_aux        additional data
   \param    los          line of sight structure
   \param    absloswfs    absorption LOS Wfs
   \param    f_mono       frequency absoprtion grid
   \param    k_grid       retrieval grid
   \param    order        polynomial order
   \param    ilow         index for lower frequency limit
   \param    iupp         index for upper frequency limit

   \author Patrick Eriksson
   \date   2000-09-15

   Adapted to MTL. 
   \date   2001-01-05
   \author Stefan Buehler

   Made it possible to have arbitrary frequency limits (ilow and iupp).
   \date   2001-01-21
   \author Patrick Eriksson
*/
void k_contabs (
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_aux,
              const LOS&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          f_mono,
              const VECTOR&          k_grid,
              const size_t&          order,
              const size_t&          ilow,
              const size_t&          iupp )
{
  // Main sizes
  const size_t  nza = los.start.size();     // number of zenith angles  
  const size_t  np  = k_grid.size();        // number of retrieval altitudes
  const size_t  npoints = order+1;          // number of off-set points
  const size_t  nv  = f_mono.size();        // number of frequencies

  // -log(p) is used as altitude variable. The minus is included to get
  // increasing values, a demand for the grid functions. 
  //  const VECTOR  lgrid=-1.0*log(k_grid);
  VECTOR  lgrid;
  p2grid( lgrid, k_grid );
  VECTOR  lplos;                      // -log of the LOS pressures

  // Indices
  // IP0 and IF0 are the index off-sets for the total K matrix
        size_t  ipoint;                    // Off-set point index
        size_t  iza;                       // Zenith angle index
        size_t  ip, ip0=0;                 // Retrieval point indices
        size_t  iv, iv0;                   // Frequency indices
        size_t  i1, iw;                    // weight indices

  // Frequency limits
  if ( ilow > iupp )
    throw logic_error("The frequency limits must be ordered.");
  if ( ilow < 0 )
    throw logic_error("The index for the lower limit must be >=0.");
  if ( iupp >= nv )
    throw logic_error(
                   "The index for the upper limit is above length of f_mono.");

  // Other variables
        VECTOR  fpoints;                   // frequencies of the off-set points
        VECTOR  b(nv);                     // fit base function
        MATRIX  is;                        // matrix for storing LOS index
        VECTOR  w;                         // weights for LOS WFs
        VECTOR  a(nv);                     // temporary vector

  // Check that the selected polynomial order
  if ( order < 0 )
    throw runtime_error("The polynomial order must be >= 0."); 

  // Set up K and additional data. Set all values of K to 0
  resize(k,nza*nv,npoints*np);
  setto( k, 0.0 );
  resize(k_names,npoints);
  resize(k_aux,npoints*np,2);

  // Calculate the frequencies of the off-set points
  if ( npoints > 1 )
    nlinspace( fpoints, f_mono[ilow], f_mono[iupp], npoints );
  else
  {
    resize( fpoints, 1 );
    fpoints[0] = ( f_mono[ilow] + f_mono[iupp] ) / 2.0;
  }  

  // The calculations
  // Loop order:
  //   1 polynomial order
  //   2 zenith angle
  //   3 retrieval altitudes 
  //   4 frequencies
  //
  // (Another loop order should be somewhat more efficient, but to keep this
  // function consistent with k_species, this loop order was selected.)
  //
  out2 << "  You have selected " << npoints << " off-set fit points.\n";
  for ( ipoint=0; ipoint<npoints; ipoint++ ) 
  {
    out2 << "  Doing point " << ipoint << "\n";

    // Fill K_NAMES and K_AUX
    {
      ostringstream os;
      os << "Continuum: " << fpoints[ipoint] << " Hz, Point " 
         << ipoint+1 << "/" << npoints;
      k_names[ipoint] = os.str();
    }
    for ( ip=0; ip<np; ip++ )
    {
       k_aux[ip0+ip][0] = k_grid[ip];
       k_aux[ip0+ip][1] = 0.0;
    }

    // Set-up base vector for the present fit point 
    setto( b, 1.0 );
    if ( npoints > 1 )
    {
      for ( ip=0; ip<npoints; ip++ )
      {
        if ( ip != ipoint )
	{ 
          for ( iv=ilow; iv<=iupp; iv++ )
            b[iv] *= (f_mono[iv]-fpoints[ip]) / ( fpoints[ipoint]-fpoints[ip]);
	}
      }
    }

    // Set frequency zenith angle index off-set to 0
    iv0 = 0;                 

    // Loop zenith angles
    out3 << "    Zenith angle nr:\n      ";
    for ( iza=0; iza<nza; iza++ ) 
    {
      if ( ((iza+1)%20)==0 )
        out3 << "\n      ";
      out3 << " " << iza; cout.flush();
      
      // Do something only if LOS has any points
      if ( los.p[iza].size() > 0 )
      {
        // Get the LOS points affected by each retrieval point
	//        lplos = -1.0 * log(los.p[iza]);
        p2grid( lplos, los.p[iza] );
        grid2grid_index( is, lplos, lgrid );

        // Loop retrieval points
        for ( ip=0; ip<np; ip++ ) 
        {
          // Check if there is anything to do
          if ( is[ip][0] >= 0 )
          {
            // Get the weights for the LOS points
	    grid2grid_weights( w, lplos, (size_t)is[ip][0], (size_t) is[ip][1],
                                                                 lgrid, ip );

            // Calculate the WFs.
            // A is di/dkappa*dkappa/dkp in a compact form.
            // This is possible as the columns of dkappa/dkp are identical.  
            setto( a, 0.0 );                     
            i1 = (size_t)is[ip][0];       // first LOS point to consider
            for ( iv=ilow; iv<=iupp; iv++ )
	    {
              for ( iw=i1; iw<=(size_t)is[ip][1]; iw++ )
                a[iv] += absloswfs[iza][iv][iw] * w[iw-i1];
              k[iv0+iv][ip0+ip] = a[iv] * b[iv];                    
	    }
          }            
        }
      }

      // Update the frequency index offset
      iv0 += nv;
    }  
    out3 << "\n";

    // Increase retrieval altitude index off-set 
    ip0 += np;
  }  
}



//// k_temp_nohydro /////////////////////////////////////////////////////////
/**
   Temparature WFs without hydrostatic eq.

   Calculates temperature 1D weighting functions WITHOUT including
   hydrostatic equilibrium. The function uses precalculated absorption
   LOS WFs, while source function LOS WFs are calculated inside the
   function.

   The expression used are described in sub-section 7 of the AUG section
   "Atmospheric WFs".

   \retval   k                 weighting function matrix
   \retval   k_names           identity name(s)
   \retval   k_aux             additional data

   \param    tag_groups        The list of tag groups
   \param    los               line of sight structure
   \param    absloswfs         absorption LOS Wfs
   \param    f_mono            frequency absorption grid
   \param    p_abs             pressure grid for abs. calculations
   \param    t_abs             temperatures at p_abs
   \param    vmrs              VMR profiles at p_abs
   \param    lines_per_tg      lines tag sorted
   \param    lineshape         lineshape specifications: function, norm, cutoff
   \param    abs               total absorption
   \param    trans             transmissions         
   \param    e_ground          ground emissivity
   \param    k_grid            retrieval grid
   \param    cont_description_names names of different continuum
                                    models
   \param    cont_description_parameters continuum parameters for the
                                         models listed in
					 cont_description_names 

   \author Patrick Eriksson
   \date   2000-09-15
*/
void k_temp_nohydro (
		           MATRIX&                        k,
		           ARRAYofstring&                 k_names,
		           MATRIX&                        k_aux,
		     const TagGroups&                     tag_groups,
		     const LOS&                           los,           
		     const ARRAYofMATRIX&                 absloswfs,
		     const VECTOR&                        f_mono,
		     const VECTOR&                        p_abs,
		     const VECTOR&                        t_abs,
		     const VECTOR&                        n2_abs,	   
		     const VECTOR&                        h2o_abs,	   
		     const ARRAYofVECTOR&                 vmrs,
		     const ARRAYofARRAYofLineRecord&      lines_per_tg,
		     const ARRAYofLineshapeSpec&          lineshape,
		     const MATRIX&                        abs,            
		     const ARRAYofMATRIX&                 trans,
		     const VECTOR&                        e_ground,
		     const VECTOR&                        k_grid,
		     const ARRAYofstring&         cont_description_names,
                     const ARRAYofVECTOR& 	  cont_description_parameters )
{
  // Main sizes
  const size_t  nza = los.start.size();     // number of zenith angles  
  const size_t  nv  = f_mono.size();        // number of frequencies
  const size_t  np  = k_grid.size();        // number of retrieval altitudes

  // -log(p) is used as altitude variable. The minus is included to get
  // increasing values, a demand for the grid functions. 
  //  const VECTOR  lgrid=-1.0*log(k_grid);
  VECTOR  lgrid;
  p2grid( lgrid, k_grid );
  VECTOR  lplos;                     // -log of the LOS pressures

  // Indices
  // IP0 and IF0 are the index off-sets for the total K matrix
        size_t  iza;                       // zenith angle index
        size_t  ip;                        // retrieval point index
        size_t  iv, iv0=0;                 // frequency indices
        size_t  i1, iw;                    // weight indices

  // Other variables
        VECTOR  t(k_grid.size());          // temperature at retrieval points
        MATRIX  abs1k;                     // absorption for t_abs+1K
        MATRIX  dabs_dt;                   // see below
 ARRAYofMATRIX  abs_dummy;                 // dummy absorption array
 ARRAYofMATRIX  sloswfs;                   // source LOS WFs
	MATRIX  is;                        // matrix for storing LOS index
        VECTOR  w;                         // weights for LOS WFs
        VECTOR  a(nv), b(nv), pl(f_mono.size());  // temporary vectors

 // The scalars are declared to be double to avoid possible numerical problems
 // when using float
        double  c,d;                       // temporary values


  // Set up K and additional data. Set all values of K to 0
  resize(k,nza*nv,np);
  setto(k, 0.0);
  resize(k_names,1);
  k_names[0] = "Temperature: no hydrostatic eq.";
  resize(k_aux,np,2);
  interpp( t, p_abs, t_abs, k_grid ); 
  for ( ip=0; ip<np; ip++ )
  {
     k_aux[ip][0] = k_grid[ip];
     k_aux[ip][1] = t[ip];
  }

  // Calculate absorption for t_abs + 1K to estimate the temperature derivative
  // dabs/dt, the temperature derivative of the absorption at k_grid
  out2 << "  Calculating absorption for t_abs + 1K\n";
  out2 << "  ----- Messages from absCalc: -----\n";
  //
  {
    VECTOR dummy(t_abs.size(),1.0);
    add(t_abs,dummy);

    absCalc( abs1k, abs_dummy, tag_groups, f_mono, p_abs, dummy, n2_abs, h2o_abs, vmrs, 
             lines_per_tg, lineshape, cont_description_names, cont_description_parameters );
  }
  resize(abs_dummy,0);
  //
  out2 << "  ----- Back from absCalc ----------\n";
  //
  add( scaled(abs,-1), abs1k );
  resize( dabs_dt, abs1k.nrows(), k_grid.size() );
  interpp( dabs_dt, p_abs, abs1k, k_grid );
  resize(abs1k,0,0);

  // Calculate source LOS WFs
  out2 << "  Calculating source LOS WFs\n";
  sourceloswfs( sloswfs, los, trans, f_mono, e_ground );

  // Determine the temperatures at the retrieval points
  out2 << "  Calculating temperature at retrieval points\n";
  interpp( t, p_abs, t_abs, k_grid );

  // The calculations
  // Loop order:
  //   1 zenith angle
  //   2 retrieval altitudes 
  //   3 frequencies
  //
  out2 << "  Calculating the weighting functions\n";
  out3 << "    Zenith angle nr:      ";
  for ( iza=0; iza<nza; iza++ ) 
  {
    if ( (iza%20)==0 )
      out3 << "\n      ";
    out3 << " " << iza; cout.flush();
    
    // Do something only if LOS has any points
    if ( los.p[iza].size() > 0 )
    {
      // Get the LOS points affected by each retrieval point
      //      lplos = -1.0 * log(los.p[iza]);
      p2grid( lplos, los.p[iza] );
      grid2grid_index( is, lplos, lgrid );

      // Loop retrieval points
      for ( ip=0; ip<np; ip++ ) 
      {
        // Check if there is anything to do
        if ( is[ip][0] >= 0 )
        {
          // Get the weights for the LOS points
	  grid2grid_weights( w, lplos, (size_t)is[ip][0], (size_t) is[ip][1], 
                                                                 lgrid, ip );

          // Calculate the WFs.
          // A is di/dsigma*dsigma/dSp in a compact form.
          // B is di/dkappa*dkappa/dkp in a compact form.
          // This is possible as the columns of dkappa/dkp are identical and  
	  // that dkappa/dkp = dsigma/dSp
          // C is just a temporary value
	  // PL is the Planck function for the present temperature value
	  //
          setto( a, 0.0 );
	  setto( b, 0.0 );                    
          c  = PLANCK_CONST / BOLTZMAN_CONST / t[ip];
          planck( pl, f_mono, t[ip] );
          i1 = (size_t)is[ip][0];       // first LOS point to consider
          //
          for ( iv=0; iv<nv; iv++ )
	  {
            for ( iw=i1; iw<=(size_t)is[ip][1]; iw++ )
	    {
              a[iv] += sloswfs[iza][iv][iw] * w[iw-i1];
              b[iv] += absloswfs[iza][iv][iw] * w[iw-i1];
	    }
            d = c * f_mono[iv];
            k[iv0+iv][ip] = a[iv] * d/t[ip] / (1-exp(-d)) * pl[iv] +
                                                      b[iv] * dabs_dt[iv][ip];
	  }
        }            
      }
    }

    // Update the frequency index offset
    iv0 += nv;
  }  
   out3 << "\n";
}




////////////////////////////////////////////////////////////////////////////
//   Workspace methods
////////////////////////////////////////////////////////////////////////////


/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Axel von Engeln
   \date   2000-01-31
*/
void wfs_tgsDefine(// WS Output:
		   TagGroups& wfs_tag_groups,
		   // Control Parameters:
		   const ARRAYofstring& tags)
{
  wfs_tag_groups = TagGroups(tags.size());

  // Each element of the array of strings tags defines one tag
  // group. Let's work through them one by one.
  for ( size_t i=0; i<tags.size(); ++i )
    {
      // There can be a comma separated list of tag definitions, so we
      // need to break the string apart at the commas.
      ARRAY<string> tag_def;

      bool go_on = true;
      string these_tags = tags[i];
      while (go_on)
	{
	  size_t n = these_tags.find(',');
	  if ( n >= these_tags.size() )
	    {
	      // There are no more commas.
	      tag_def.push_back(these_tags);
	      go_on = false;
	    }
	  else
	    {
	      tag_def.push_back(these_tags.substr(0,n));
	      these_tags.erase(0,n+1);
	    }
	}

      // Unfortunately, MTL conatains a bug that leads to all elements of
      // the outer ARRAY of an ARRAY<ARRAY>> pointing to the same data
      // after creation. So we need to fix this explicitly:
      wfs_tag_groups[i] = ARRAY<OneTag>();

      for ( size_t s=0; s<tag_def.size(); ++s )
	{
	  // Remove leading whitespace, if there is any:
	  while ( ' '  == tag_def[s][0] ||
		  '\t' == tag_def[s][0]    )	tag_def[s].erase(0,1);

	  OneTag this_tag(tag_def[s]);

	  // Safety check: For s>0 check that the tags belong 
          // to the same species.
	  if (s>0)
	    if ( wfs_tag_groups[i][0].Species() != this_tag.Species() )
	      throw runtime_error(
                       "Tags in a tag group must belong to the same species.");

	  wfs_tag_groups[i].push_back(this_tag);
	}
    }

  // Print list of tag groups to the most verbose output stream:
  out3 << "  Defined weighting function tag groups:";
  for ( size_t i=0; i<wfs_tag_groups.size(); ++i )
    {
      out3 << "\n  " << i << ":";
      for ( size_t s=0; s<wfs_tag_groups[i].size(); ++s )
	{
	  out3 << " " << wfs_tag_groups[i][s].Name();
	}
    }
  out3 << '\n';
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-03-30
*/
void absloswfsCalc (
                    ARRAYofMATRIX&   absloswfs,
	      const int&             emission,
              const LOS&             los,   
              const ARRAYofMATRIX&   source,
              const ARRAYofMATRIX&   trans,
              const VECTOR&          y,
              const VECTOR&          y_space,
              const VECTOR&          f_mono,
              const VECTOR&          e_ground,
              const Numeric&         t_ground )
{
  if ( emission == 0 )
    absloswfs_tau ( absloswfs, los, f_mono );
  else
    absloswfs_rte ( absloswfs, los, source, trans, y, y_space, f_mono,     
                                                          e_ground, t_ground );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-03-30
*/
void absloswfsTau (
                    ARRAYofMATRIX&   absloswfs,
	      const int&             emission,
	      const LOS&             los,
              const VECTOR&          f_mono )
{
  if ( emission != 0 )
    throw runtime_error(
      "The function yTau can only be used when emission is neglected.");

  absloswfs_tau ( absloswfs, los, f_mono );

}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kSpecies (
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_aux,
              const LOS&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          p_abs,
              const VECTOR&          t_abs,             
              const TagGroups&       wfs_tgs,
              const ARRAYofMATRIX&   abs_per_tg,
              const ARRAYofVECTOR&   vmrs,
              const VECTOR&          k_grid,
              const string&          tag,
              const string&          unit )
{
  if ( wfs_tgs.size() != abs_per_tg.size() )
    throw runtime_error( "Lengths of wfs_tgs and abs_per_tg do not match." ); 

  ARRAYofstring  tag_name(1);
  tag_name[0] = tag;

  ARRAYofsizet   tg_nr; 
  get_tagindex_for_strings( tg_nr, wfs_tgs, tag_name );
  
  k_species( k, k_names, k_aux, los, absloswfs, p_abs, t_abs, 
                              wfs_tgs, abs_per_tg, vmrs, k_grid, tg_nr, unit );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kSpeciesAll (
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_aux,
              const LOS&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          p_abs,
              const VECTOR&          t_abs,             
              const TagGroups&       wfs_tgs,
              const ARRAYofMATRIX&   abs_per_tg,
              const ARRAYofVECTOR&   vmrs,
              const VECTOR&          k_grid,
              const string&          unit )
{
  const size_t  ntg = wfs_tgs.size();     // number of retrieval tags
  ARRAYofsizet  tg_nr(ntg);

  if ( ntg != abs_per_tg.size() )
    throw runtime_error( "Lengths of wfs_tgs and abs_per_tg do not match." ); 
  
  for ( size_t i=0; i<ntg; i++ )
    tg_nr[i] = i;

  k_species( k, k_names, k_aux, los, absloswfs, p_abs, t_abs, 
                              wfs_tgs, abs_per_tg, vmrs, k_grid, tg_nr, unit );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kContAbs (
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_aux,
              const LOS&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          f_mono,
              const VECTOR&          k_grid,
              const int&             order )
{
  k_contabs( k, k_names, k_aux, los, absloswfs, f_mono, k_grid, order,
                                                          0, f_mono.size()-1 );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-01-21
*/
void kContAbsSpecifiedLimits (
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_aux,
              const LOS&             los,           
              const ARRAYofMATRIX&   absloswfs,
              const VECTOR&          f_mono,
              const VECTOR&          k_grid,
		    const int&             order,
              const Numeric&         f_low,
              const Numeric&         f_high )
{
  if ( f_low > f_high )
    throw runtime_error(
       "The lower frequency limit has a higher value than the upper limit." ); 

  const INDEX   nf = f_mono.size();
        INDEX   i_low, i_high;
  for( i_low=0; i_low<nf && f_mono[i_low] < f_low; i_low++ )
    {}
  if ( i_low == nf )
    throw runtime_error(
       "The lower frequency limit is above all values of f_mono." ); 

  for( i_high=i_low; i_high<nf && f_mono[i_high] <= f_high; i_high++ )
    {}
  if ( i_high == 0 )
    throw runtime_error(
       "The upper frequency limit is below all values of f_mono." ); 
  i_high--;

  k_contabs( k, k_names, k_aux, los, absloswfs, f_mono, k_grid, order,
                                                              i_low, i_high );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kTempNoHydro (
                    MATRIX&                         k,
                    ARRAYofstring&                  k_names,
                    MATRIX&                         k_aux,
                    const TagGroups&                tag_groups,
		    const LOS&                      los,           
		    const ARRAYofMATRIX&            absloswfs,
		    const VECTOR&                   f_mono,
		    const VECTOR&                   p_abs,
		    const VECTOR&                   t_abs,
		    const VECTOR&                   n2_abs,
		    const VECTOR&                   h2o_abs,
		    const ARRAYofVECTOR&            vmrs,
		    const ARRAYofARRAYofLineRecord& lines_per_tg,
		    const ARRAYofLineshapeSpec&     lineshape,
		    const MATRIX&                   abs,            
		    const ARRAYofMATRIX&            trans,
		    const VECTOR&                   e_ground,
		    const VECTOR&                   k_grid,
		    const ARRAYofstring&            cont_description_names,
		    const ARRAYofVECTOR& 	    cont_description_parameters )
{
  k_temp_nohydro( k, k_names, k_aux, tag_groups, los, absloswfs, f_mono, p_abs,
                  t_abs, n2_abs, h2o_abs, vmrs, lines_per_tg, lineshape, 
                  abs, trans, e_ground, k_grid,
		  cont_description_names, cont_description_parameters );  
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?

   Adapted to MTL.
   \date   2001-01-06
   \author Stefan Buehler
 */
void kManual(
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_aux,
              const VECTOR&          y0,      
              const VECTOR&          y,
              const string&          name,
              const Numeric&         delta,
              const Numeric&         grid,
              const Numeric&         apriori )
{
  assert( y.size()==y0.size() );
  VECTOR dummy(y.size());		// To store intermediate result.

  // Put y-y0 in dummy:
  add( y, scaled(y0,-1), dummy );

  // Make k one-column matrix of the right size:
  resize( k, y.size(), 1 );

  // Divide dummy by delta and copy it to first column of k:
  copy( scaled(dummy,1./delta), columns(k)[0]);
  
  resize(k_names,1);
  k_names[0] = name;
  resize(k_aux,1,2);
  k_aux[0][0] = grid;
  k_aux[0][1] = apriori;
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?

   Adapted to MTL.
   \date   2001-01-06
   \author Stefan Buehler
 */
void kDiffHSmall(
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_aux,
              const Hmatrix&         h1,      
              const Hmatrix&         h2,      
              const VECTOR&          y,
              const string&          name,
              const Numeric&         delta,
              const Numeric&         grid,
              const Numeric&         apriori )
{
  VECTOR y1, y2;
  h_apply( y1, h1, y );
  h_apply( y2, h2, y );

  assert( y1.size()==y2.size() );
  VECTOR dummy(y1.size());

  // Put y2-y1 in dummy:
  add( y2, scaled(y1,-1), dummy );

  // Make k one-column matrix of the right size:
  resize( k, y1.size(), 1 );
  // Divide dummy by delta and copy it to first column of k:
  copy( scaled(dummy,1./delta), columns(k)[0]);

  resize(k_names,1);
  k_names[0] = name;
  resize(k_aux,1,2);
  k_aux[0][0] = grid;
  k_aux[0][1] = apriori;
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?

   Adapted to MTL.
   \date   2001-01-06
   \author Stefan Buehler
 */
void kDiffHFast(
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_aux,
              const Hmatrix&         h1,      
              const Hmatrix&         h2,      
              const VECTOR&          y,
              const string&          name,
              const Numeric&         delta,
              const Numeric&         grid,
              const Numeric&         apriori )
{
  VECTOR yd;
  Hmatrix hd;
  h_diff( hd, h2, h1 );
  h_apply( yd, hd, y );
  
  // Make k one-column matrix of the right size:
  resize( k, yd.size(), 1 );

  // Divide yd by delta and copy it to first column of k:
  copy( scaled(yd,1./delta), columns(k)[0]);

  resize(k_names,1);
  k_names[0] = name;
  resize(k_aux,1,2);
  k_aux[0][0] = grid;
  k_aux[0][1] = apriori;
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-21
 */
void kPointingOffSet(
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_aux,
              const Numeric&         z_plat,
              const VECTOR&          za_pencil,
              const Numeric&         l_step,
              const VECTOR&          p_abs,
              const VECTOR&          z_abs,
              const VECTOR&          t_abs,
              const VECTOR&          f_mono,
              const int&             refr,
              const int&             refr_lfac,
              const VECTOR&          refr_index,
              const Numeric&         z_ground,
              const Numeric&         r_geoid,
              const MATRIX&          abs,
  	      const int&             emission,
              const VECTOR&          y_space,
              const VECTOR&          e_ground,
              const Numeric&         t_ground,
              const VECTOR&          y,
              const Numeric&         delta )
{
  // Create new zenith angle grid
  const INDEX  nza = za_pencil.size();
  VECTOR za_new( nza );
  copy( za_pencil, za_new );
  add( za_new, VECTOR(nza,delta), za_new );

  out2 << "  ----- Messages from losCalc: --------\n";
  LOS    los;
  VECTOR z_tan;
  losCalc( los, z_tan, z_plat, za_new, l_step, p_abs, z_abs, refr, refr_lfac, 
                                              refr_index, z_ground, r_geoid );
  out2 << "  -------------------------------------\n";

  ARRAYofMATRIX source, trans;
  out2 << "  ----- Messages from sourceCalc: -----\n";
  sourceCalc( source, emission, los, p_abs, t_abs, f_mono );
  out2 << "  -------------------------------------\n";
  out2 << "  ----- Messages from transCalc: ------\n";
  transCalc( trans, los, p_abs, abs );
  out2 << "  -------------------------------------\n";
  VECTOR y_new;
  out2 << "  ----- Messages from yRte: -----------\n";
  yCalc( y_new, emission, los, f_mono, y_space, source, trans, 
                                                         e_ground, t_ground );
  out2 << "  -------------------------------------\n";

  // Make k one-column matrix of the right size:
  const INDEX   nv = y.size();
  resize( k, nv, 1 );

  // k = (y_new - y) / delta
  VECTOR dummy( nv );
  add( y_new, scaled(y,-1), dummy );
  copy( scaled(dummy,1/delta), columns(k)[0] );

  resize( k_names, 1 );
  k_names[0] = "Pointing: off-set";
  resize( k_aux, 1, 2 );
  setto( k_aux, 0.0 );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-21
 */
void kCalibration(
                    MATRIX&          k,
                    ARRAYofstring&   k_names,
                    MATRIX&          k_aux,
              const VECTOR&          y,
              const VECTOR&          y_ref,
       	      const string&          name )
{
  const INDEX   ny = y.size();
  const INDEX   nf = y_ref.size();
  const INDEX   nza = INDEX( floor(ny/nf) );

  if ( nza*nf != ny )
    throw runtime_error("The length of y_ref does not match the length of y");
  
  // Make k one-column matrix of the right size:
  resize( k, ny, 1 );

  // k = y - y_ref
  INDEX j,i,i0;
  for ( j=0; j<nza; j++ )    
  {
    i0 = j*nf;
    for ( i=0; i<nf; i++ )
      k[i0+i][0] = y[i0+i] - y_ref[i];
  }

  resize( k_names, 1 );
  k_names[0] = "Calibration: scaling";
  resize( k_aux, 1, 2 );
  setto( k_aux, 0.0 );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kxInit (
                    MATRIX&          kx,
                    ARRAYofstring&   kx_names,
                    ARRAYofsizet&    kx_lengths,
                    MATRIX&          kx_aux )
{
  resize( kx,         0, 0 );
  resize( kx_names,   0    );
  resize( kx_lengths, 0    );
  resize( kx_aux,     0, 0 );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kbInit (
                    MATRIX&          kb,
                    ARRAYofstring&   kb_names,
                    ARRAYofsizet&    kb_lengths,
                    MATRIX&          kb_aux )
{
  kxInit( kb, kb_names, kb_lengths, kb_aux );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kxAppend (
		    MATRIX&          kx,
		    ARRAYofstring&   kx_names,
		    ARRAYofsizet&    kx_lengths,
		    MATRIX&          kx_aux,
              const MATRIX&          k,
              const ARRAYofstring&   k_names,
              const MATRIX&          k_aux )
{
  k_append( kx, kx_names, kx_lengths, kx_aux, k, k_names, k_aux );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kbAppend (
		    MATRIX&          kb,
		    ARRAYofstring&   kb_names,
		    ARRAYofsizet&    kb_lengths,
		    MATRIX&          kb_aux,
              const MATRIX&          k,
              const ARRAYofstring&   k_names,
              const MATRIX&          k_aux )
{
  k_append( kb, kb_names, kb_lengths, kb_aux, k, k_names, k_aux );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-21
*/
void kxAllocate (
                    MATRIX&          kx,
                    ARRAYofstring&   kx_names,
                    ARRAYofsizet&    kx_lengths,
                    MATRIX&          kx_aux,
              const VECTOR&          y,
              const string&          y_name,
              const int&             ni,
              const int&             nx )
{
  const INDEX  ny = y.size();

  out2 << "  Allocates memory for a weighting function matrix having:\n" 
       << "    " << ny << " frequency points\n"
       << "    " << nx << " state variables\n"
       << "    " << ni << " retrieval identities\n";

  resize( kx,         ny, nx );
  resize( kx_names,   ni     );
  resize( kx_lengths, ni     );
  resize( kx_aux,     nx, 2  );

  for ( int i=0; i<ni; i++ )
    kx_lengths[i] = 0;
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-21
*/
void kbAllocate (
                    MATRIX&          kb,
                    ARRAYofstring&   kb_names,
                    ARRAYofsizet&    kb_lengths,
                    MATRIX&          kb_aux,
              const VECTOR&          y,
              const string&          y_name,
              const int&             ni,
              const int&             nx )
{
  kxAllocate( kb, kb_names, kb_lengths, kb_aux, y, y_name, ni, nx );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-21
*/
void kxPutInK (
		    MATRIX&          kx,
		    ARRAYofstring&   kx_names,
		    ARRAYofsizet&    kx_lengths,
		    MATRIX&          kx_aux,
              const MATRIX&          k,
              const ARRAYofstring&   k_names,
              const MATRIX&          k_aux )
{
  const INDEX  ny  = kx.nrows();
  const INDEX  nx  = kx.ncols();
  const INDEX  ni  = kx_names.size();
  const INDEX  nx2 = k.ncols();
  const INDEX  ni2 = k_names.size();

  if ( ny != k.nrows() )
    throw runtime_error("The two WF matrices have different number of rows.");

  // Determine old length of x and number of identities
  INDEX  ni1, nx1=0;
  for( ni1=0; ni1<ni && kx_lengths[ni1]>0; ni1++ )
    nx1 += kx_lengths[ni1];
  if ( ni1 == ni )
    throw runtime_error("All retrieval/error identity positions used.");

  // Check if new data fit
  if ( (ni1+ni2) > ni )
  {
    ostringstream os;
    os << "There is not room for so many retrieval/error identities.\n" 
       << (ni1+ni2)-ni << " positions are missing";
      throw runtime_error(os.str());
  }
  if ( (nx1+nx2) > nx )
  {
    ostringstream os;
    os << "The k-matrix does not fit in kx.\n" 
       << (nx1+nx2)-nx << " columns are missing";
      throw runtime_error(os.str());
  }

  out2 << "  Putting k in kx as\n" 
       << "    identity " << ni1 << " - " << ni1+ni2-1 << "\n"
       << "      column " << nx1 << " - " << nx1+nx2-1 << "\n";
 
  // Calculate the vector length for each identity in K
  size_t l = (size_t) floor(nx2/ni2);

  // Move K to Kx
  copy( k,       kx.sub_matrix( 0, ny, nx1, nx1+nx2 ) );
  copy( k_aux,   kx_aux.sub_matrix( nx1, nx1+nx2, 0, 2 ) );
  for ( INDEX iri=0; iri<ni2; iri++ )
  {
    kx_names[ni1+iri]   = k_names[iri];
    kx_lengths[ni1+iri] = l;
  }   
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-21
*/
void kbPutInK (
		    MATRIX&          kb,
		    ARRAYofstring&   kb_names,
		    ARRAYofsizet&    kb_lengths,
		    MATRIX&          kb_aux,
              const MATRIX&          k,
              const ARRAYofstring&   k_names,
              const MATRIX&          k_aux )
{
  kxPutInK( kb, kb_names, kb_lengths, kb_aux, k, k_names, k_aux );
}
