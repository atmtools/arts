/* Copyright (C) 2000, 2001 Patrick Eriksson <patrick@rss.chalmers.se>

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

#include <math.h>
#include "arts.h"
#include "auto_md.h"
#include "matpackI.h"
#include "messages.h"          
#include "auto_wsv.h"          
#include "math_funcs.h"          
#include "atm_funcs.h"          
#include "los.h"
extern const Numeric PLANCK_CONST;
extern const Numeric BOLTZMAN_CONST;



////////////////////////////////////////////////////////////////////////////
//   Function(s) to join two WF matrices and related variables
////////////////////////////////////////////////////////////////////////////

//// k_append ////////////////////////////////////////////////////////////////
/**
   Appends the K matrix to either Kx or Kb. 

   Output matrix kx is resized!

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
		    Matrix&          kx,
		    ArrayOfString&   kx_names,
		    ArrayOfIndex&    kx_lengths,
		    Matrix&          kx_aux,
              const Matrix&          k,
              const ArrayOfString&   k_names,
              const Matrix&          k_aux )
{
  // Size of Kx and K
  const Index  ny1  = kx.nrows();         // length of measurement vector (y)
  const Index  nx1  = kx.ncols();         // length of state vector (x)
  const Index  nri1 = kx_names.nelem();    // number of retrieval identities
  const Index  ny2  = k.nrows();  
  const Index  nx2  = k.ncols();  
  const Index  nri2 = k_names.nelem();
  Index        iri;

  // Make copy of Kx data.
  Matrix          ktemp(kx);
  ArrayOfString   ktemp_names(kx_names);
  ArrayOfIndex    ktemp_lengths(kx_lengths);
  Matrix          ktemp_aux(kx_aux);
  
//   cout << "ktemp =\n" << ktemp << "\n\n\n";
//   cout << "ktemp_names =\n" << ktemp_names << "\n\n\n";
//   cout << "ktemp_lengths =\n" << ktemp_lengths << "\n\n\n";
//   cout << "ktemp_aux =\n" << ktemp_aux << "\n\n\n";

  if ( nx1 > 0 )
  {
    if ( ny1 != ny2 )
      throw runtime_error(
            "The two WF matrices have different number of rows." ); 
  }

  // Reallocate the Kx data
  kx.resize(         ny2,       nx1+nx2 );
  kx_names.resize(   nri1+nri2          );
  kx_lengths.resize( nri1+nri2          );
  kx_aux.resize(     nx1+nx2,   2       );

  // Move Ktemp to Kx
  if ( nx1 > 0 )
  {
    kx( Range(joker), Range(0,nx1) )     = ktemp;
    kx_aux( Range(0,nx1), Range(joker) ) = ktemp_aux;
    // For the Array types kx_names and kx_lenghts we cannot use Range.  
    for ( iri=0; iri<nri1; iri++ )
    {
      kx_lengths[iri]  = ktemp_lengths[iri];
      kx_names[iri]    = ktemp_names[iri];
    }
  }

  // Calculate the vector length for each identity in K
  Index l = (Index) floor(nx2/nri2);

  // Move K to Kx
  kx( Range(joker), Range(nx1,nx2) )     = k;
  kx_aux( Range(nx1,nx2), Range(joker) ) = k_aux;
  // For the Array types kx_names and kx_lenghts we cannot use Range.  
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
	      Vector&   grid,
        const Vector&   pgrid )
{
  grid.resize( pgrid.nelem() );
  grid = pgrid;			// Matpack can copy the contents of
				// vectors like this. The dimensions
				// must be the same! 
  transform( grid, log, grid );	// Matpack has the same transform
				// functionality like the old
				// vecmat. THE ORDER OF THE ARGUMENTS
				// IS DIFFERENT, THOUGH (output first).
  grid *= -1;			// Mulitply all elements with -1.
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
                    Matrix&   is,
              const Vector&   x0,
              const Vector&   xp )
{
  const Index n0 = x0.nelem();        // length of original grid
  const Index np = xp.nelem();        // length if retrieval grid
  Index       i0, ip;                // counter for each grid

  // Resize is and set all values to -1
  is.resize( np, 2 );
  is = -1.0;			// Matpack can set all elements like this.
 
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
      is(0,0) = (Numeric) i0;
      for ( ; (i0<n0) && (x0[i0]<xp[1]); i0++ ) {}
      is(0,1) = (Numeric) (i0 - 1);
    }

    // Points inside XP
    for ( ip=1; ip<(np-1); ip++ )
    {
      i0 = 0;
      for ( ; (i0<n0) && (x0[i0]<=xp[ip-1]); i0++ ) {}
      if ( (i0<n0) && (x0[i0]<xp[ip+1]) )
      {
        is(ip,0) = (Numeric) i0;
        for ( ; (i0<n0) && (x0[i0]<xp[ip+1]); i0++ ) {}
        is(ip,1) = (Numeric) (i0 - 1);
      }
    }

    // Last point of XP
    i0 = 0;
    for ( ; (i0<n0) && (x0[i0]<=xp[np-2]); i0++ ) {}
    if ( (i0<n0) && (x0[i0]<=xp[np-1]) )
    {
      is(np-1,0) = (Numeric) i0;
      for ( ; (i0<n0) && (x0[i0]<=xp[np-1]); i0++ ) {}
      is(np-1,1) = (Numeric) (i0 - 1);
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
                    Vector&   w,
              const Vector&   x0,
              const Index&   i1,
              const Index&   i2,
              const Vector&   xp,
              const Index&   ip )
{
  const Index   np = xp.nelem();        // number of retrieval points
  const Index   nw = i2-i1+1;          // number of LOS points affected
  Index         i;

  // Reallocate w
  w.resize(nw);

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
                    Matrix&   k,
                    Vector    y,
              const Index&   start_index,
	      const Index&   stop_index,   // this variable is used by
              const Numeric&  lstep,        // absloswfs_down
              const Matrix&   tr,
              const Matrix&   s,
              const Index&    ground,
              const Vector&   e_ground,
              const Vector&   y_ground )
{
  const Index   nf = tr.nrows();     // number of frequencies
  Index         iv;                  // frequency index
  Index         q;                   // LOS index (same notation as in AUG)
  Vector        t(nf);               // transmission to the sensor
  
  t = 1.0;		// Set all elements of t to 1.

  if ( ground && ((ground-1==stop_index) || (ground-1==start_index)) )
    throw logic_error("The ground cannot be at one of the end points."); 

  // Resize K
  k.resize( nf, start_index+1 );

  // We start here at the LOS point closest to the sensor, that is,
  // reversed order compared to RTE_ITERATE  

  // The LOS point closest to the sensor 
  q  = stop_index;
  for ( iv=0; iv<nf; iv++ )    
  {
    t[iv]   *= tr(iv,q);
    y[iv]    = (y[iv]-s(iv,q)*(1.0-tr(iv,q)))/tr(iv,q);
    k(iv,q) = (-lstep/2)*(y[iv]-s(iv,q))*t[iv];
  }

  // Points inside the LOS
  for ( q=stop_index+1; q<start_index; q++ )
  {
    // Not a ground point
    if ( q != (ground-1) )
    {
      for ( iv=0; iv<nf; iv++ )    
      {
        y[iv]    = (y[iv]-s(iv,q)*(1.0-tr(iv,q)))/tr(iv,q);
        k(iv,q) = (-lstep/2)*( 2*(y[iv]-s(iv,q))*tr(iv,q) + 
                                             s(iv,q) - s(iv,q-1) ) *t[iv];
        t[iv]   *= tr(iv,q); 
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
                            s(iv,q)*(1.0-tr(iv,q))*(1-e_ground[iv]) ) / 
                                             ( tr(iv,q) * (1-e_ground[iv]) );
        k(iv,q) = (-lstep/2)*( 2*(y[iv]-s(iv,q))*tr(iv,q)*(1-e_ground[iv])+ 
                                s(iv,q)*(1-e_ground[iv]) + 
                             e_ground[iv]*y_ground[iv] - s(iv,q-1) ) * t[iv];
        t[iv]   *= tr(iv,q) * (1-e_ground[iv]); 
      }
    }
  }

  // The LOS point furthest away from the sensor
  q = start_index;
  for ( iv=0; iv<nf; iv++ )    
    k(iv,q)  = (-lstep/2)*(y[iv]-s(iv,q-1))*t[iv];

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
                    Matrix&   k,
                    Vector    y,              // = y_q^q
              const Vector&   y_space,             
              const Index&   start_index,
              const Numeric&  lstep,
              const Matrix&   tr,
              const Matrix&   s,
              const Index&    ground,
              const Vector&   e_ground )
{
  const Index nf = tr.nrows();     // number of frequencies
  Index       iv;                  // frequency index
  Vector      t1q;                 // transmission tangent point - q squared
  Vector      tqn(nf,1.0);         // transmission q - sensor
  Index       q;                   // LOS index (same notation as in AUG)
  Numeric     tv, tv1;             // transmission value for q and q-1
  Vector      yn(y_space);         // = y_space

  // Matpack can initialize a new Vector from another Vector. See how
  // yn is initialized from y_space.

  // Resize K
  k.resize( nf, start_index+1 );

  // Calculate the total transmission
  bl( t1q, start_index, start_index, tr, ground, e_ground );

  // We start at the outermost point
  q  = start_index;       
  for ( iv=0; iv<nf; iv++ )    
  {
    tv1      = tr(iv,q-1);
    t1q[iv] /= tv1*tv1;
    tqn[iv] *= tv1;
    y[iv]    = ( y[iv] - s(iv,q-1)*(1-tv1)*(1+t1q[iv]*tv1) ) / tv1;
    k(iv,q)  = (-lstep/2)*( ( 2*yn[iv] + s(iv,q-1)*(1-2*tv1) ) * 
                              t1q[iv]*tv1 + y[iv] - s(iv,q-1) )*tv1;
  }

  // Points inside the LOS
  for ( q=start_index-1; q>0; q-- )
  {
    for ( iv=0; iv<nf; iv++ )    
    {
      tv1      = tr(iv,q-1);    
      tv       = tr(iv,q);
      t1q[iv] /= tv1*tv1;
      y[iv]    = ( y[iv] - s(iv,q-1)*(1-tv1)*(1+t1q[iv]*tv1) ) / tv1;
      k(iv,q) = (-lstep/2) * ( 
           ( 4*(yn[iv]-s(iv,q))*tv1*tv + 3*(s(iv,q)-s(iv,q-1))*tv1 + 
                                        2*s(iv,q-1) ) * t1q[iv]*tv1 + 
             2*(y[iv]-s(iv,q-1))*tv1 + s(iv,q-1) - s(iv,q) ) * tqn[iv];
      tqn[iv] *= tv1;
      yn[iv]   = yn[iv]*tv + s(iv,q)*(1-tv);
    } 
  }

  // The tangent or ground point
  for ( iv=0; iv<nf; iv++ )    
    k(iv,0)  = (-lstep/2)*( (2*yn[iv]*tv1+s(iv,0)*(1-2*tv1))*t1q[iv] +
                       y[iv] - s(iv,0) ) * tqn[iv];

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
                    Matrix&   k,
              const Vector&   y,
              const Vector&   y_space,
              const Index&   start_index,
              const Index&   stop_index,
              const Numeric&  lstep,
              const Matrix&   tr,
              const Matrix&   s,
              const Index&    ground,
              const Vector&   e_ground,
              const Vector&   y_ground )
{
  const Index   nf = tr.nrows(); // number of frequencies
        Index   iv;              // frequency index
        Index   q;               // LOS index (same notation as in AUG)
        Vector   y0(nf);          // see below
        Matrix   k2;              // matrix for calling other LOS WFs functions
        Vector   tr0;             // see below

  // Resize K
  k.resize( nf, start_index+1 );

  // Calculate Y0, the intensity reaching the platform altitude at the far
  // end of LOS, that is, from above.
  y0 = y_space;			// Matpack can copy the contents of
				// vectors like this. The dimensions
				// must be the same! 
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
      k(iv,q) = k2(iv,q);
  }

  // The indeces above STOP_INDEX are handled by the 1pass function.
  // The transmission below STOP_INDEX must here be considered.
  absloswfs_rte_1pass( k2, y0, start_index, stop_index, lstep, tr, s, 
                     ground, e_ground, y_ground );
  for ( iv=0; iv<nf; iv++ )
  {
    for ( q=stop_index+1; q<=start_index; q++ )
      k(iv,q) = k2(iv,q)*tr0[iv];
  }

  // The platform altitude must be treated seperately
  //
  // Calculate the intensity generated below point q-1, YQQ
  Vector yqq(nf);
  rte( yqq, stop_index-1, stop_index-1, tr, s, Vector(nf,0.0), 
                                            ground, e_ground, y_ground );
  //
  // Y0 is moved one step upwards and TR0 one step downwards
  q = stop_index; 
  for ( iv=0; iv<nf; iv++ )
  {
    tr0[iv] /= tr(iv,q-1)*tr(iv,q-1);    
    y0[iv]   = ( y0[iv] - s(iv,q)*(1-tr(iv,q)) ) / tr(iv,q);
    k(iv,q) = (-lstep/2)*( (3*(y0[iv]-s(iv,q))*tr(iv,q-1)*tr(iv,q) + 
                  2*(s(iv,q)-s(iv,q-1))*tr(iv,q-1) + s(iv,q-1) )*tr0[iv] + 
                                         yqq[iv] - s(iv,q-1) )*tr(iv,q-1);
  }
}



//// absloswfs_rte ///////////////////////////////////////////////////////////
/**
   \author Patrick Eriksson
   \date   2000-09-15
*/
void absloswfs_rte (
                    ArrayOfMatrix&   absloswfs,
              const Los&             los,   
              const ArrayOfMatrix&   source,
              const ArrayOfMatrix&   trans,
              const Vector&          y,
              const Vector&          y_space,
              const Vector&          f_mono,
              const Vector&          e_ground,
              const Numeric&         t_ground )
{
  const Index  nza = los.start.nelem();   // number of zenith angles  
  const Index  nf  = f_mono.nelem();      // number of frequencies  
        Vector  yp(nf);                   // part of Y
        Index  iy0=0;                    // y index

  // Set up vector for ground blackbody radiation
  Vector   y_ground(f_mono.nelem()); 
  if ( any_ground(los.ground) )
    planck( y_ground, f_mono, t_ground );

  // Resize the LOS WFs array
  absloswfs.resize( nza );

  // Loop zenith angles
  out3 << "    Zenith angle nr:\n      ";
  for (Index i=0; i<nza; i++ ) 
  {
    if ( ((i+1)%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    // Do something only if LOS has any points
    if ( los.p[i].nelem() > 0 )
    {

      // Pick out the part of Y corresponding to the present zenith angle
      yp = y[Range(iy0,nf)];

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
                    ArrayOfMatrix&   absloswfs,
	      const Los&             los,
	      const Vector&          f_mono )
{
  const Index  nza = los.start.nelem();   // number of zenith angles  
  const Index  nf  = f_mono.nelem();      // number of frequencies  
  Numeric      kw, kw2;
  Index        row, col, np;

  // Resize the LOS WFs array
  absloswfs.resize( nza );

  // Loop zenith angles
  out3 << "    Zenith angle nr:\n      ";
  for (Index i=0; i<nza; i++ ) 
  {
    if ( ((i+1)%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();

    np = los.start[i] + 1;

    absloswfs[i].resize( nf, np );
    
    // Do something only if LOS has any points
    if ( los.p[i].nelem() > 0 )
    {

      // Upward looking (=single pass)
      if ( los.stop[i]==0 )
      {
        kw  = los.l_step[i] / 2.0;
        kw2 = los.l_step[i];
        for( row=0; row<nf; row++ )
	{
          absloswfs[i](row,0) = kw;
          absloswfs[i](row,np-1) = kw;
          for( col=1; col<np-1; col++ )
            absloswfs[i](row,col) = kw2;
        }
      }

      // 1D limb sounding
      else if ( los.start[i] == los.stop[i] )
      {
        kw  = los.l_step[i];
        kw2 = los.l_step[i] * 2.0;
        for( row=0; row<nf; row++ )
	{
          absloswfs[i](row,0) = kw;
          absloswfs[i](row,np-1) = kw;
          for( col=1; col<np-1; col++ )
            absloswfs[i](row,col) = kw2;
        }
      }

      // 1D downward looking
      else 
      {
        kw  = los.l_step[i];
        kw2 = los.l_step[i] * 2.0;
        for( row=0; row<nf; row++ )
	{
          absloswfs[i](row,0) = kw;
          absloswfs[i](row,np-1) = kw / 2.0;
          absloswfs[i](row,los.stop[i]) = kw * 1.5;
          for( col=1; col<los.stop[i]; col++ )
            absloswfs[i](row,col) = kw2;
          for( col=los.stop[i]+1; col<np-1; col++ )
            absloswfs[i](row,col) = kw;
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
                    Matrix&   k,
              const Index&   start_index,
	      const Index&   stop_index,   // this variable is used by 1D down
              const Matrix&   tr,
              const Index&    ground,
              const Vector&   e_ground )
{
  const Index   nf = tr.nrows();     // number of frequencies
        Index   iv;                  // frequency index
        Vector   t(nf,1.0);           // transmission to the sensor
        Index   q;                   // LOS index (same notation as in AUG) 

  if ( ground && ((ground-1==stop_index) || (ground-1==start_index)) )
    throw logic_error("The ground cannot be at one of the end points."); 

  // Resize K
  k.resize( nf, start_index+1 );

  // We start here at the LOS point closest to the sensor, that is,
  // reversed order compared to RTE_ITERATE  

  // The LOS point closest to the sensor 
  q  = stop_index;
  for ( iv=0; iv<nf; iv++ )    
    k(iv,q)  = ( 1 - tr(iv,q) ) / 2;

  // Points inside the LOS
  for ( q=stop_index+1; q<start_index; q++ )
  {
    // Not a ground point
    if ( q != ground )
    {
      for ( iv=0; iv<nf; iv++ )    
      {
        k(iv,q) = ( 1 - tr(iv,q-1)*tr(iv,q) ) * t[iv] / 2;
        t[iv]  *= tr(iv,q); 
      }
    }
    // A ground point
    else
    {
      out1 << "WARNING: The function sourceloswfs_1pass not tested "
              "for ground reflections\n";

      for ( iv=0; iv<nf; iv++ )    
      {
        k(iv,q) = ( (1-tr(iv,q))*(1-e_ground[iv])*tr(iv,q-1) + 1 - 
                                                     tr(iv,q-1) ) * t[iv] / 2;
        t[iv]  *= tr(iv,q)*(1-e_ground[iv]); 
      }
    }
  }

  // The LOS point furthest away from the sensor
  q = start_index;
  for ( iv=0; iv<nf; iv++ )    
    k(iv,q)  = ( 1 - tr(iv,q-1) ) * t[iv] / 2;
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
                    Matrix&   k,
              const Index&   start_index,
              const Matrix&   tr,
              const Index&    ground,
              const Vector&   e_ground )
{
  const Index   nf = tr.nrows();      // number of frequencies
        Index   iv;                  // frequency index
        Vector   t1q;                 // transmission tangent point - q squared
        Vector   tqn(nf,1);           // transmission q - sensor
        Index   q;                   // LOS index (same notation as in AUG)

  // Calculate the total transmission
  bl( t1q, start_index, start_index, tr, ground, e_ground );

  // Resize K
  k.resize( nf, start_index+1 );

  // We start at the outermost point
  q  = start_index;       
  for ( iv=0; iv<nf; iv++ )    
  {
    t1q[iv] /= tr(iv,q-1) * tr(iv,q-1);
    k(iv,q)  = ( (1-tr(iv,q-1))*t1q[iv]*tr(iv,q-1) + 1 - tr(iv,q-1) )/2;
  }

  // Points inside the LOS
  for ( q=start_index-1; q>0; q-- )
  {
    for ( iv=0; iv<nf; iv++ )    
    {
      t1q[iv]  /= tr(iv,q-1) * tr(iv,q-1);
      k(iv,q)  = ( (1-tr(iv,q-1)*tr(iv,q))*t1q[iv]*tr(iv,q-1)*
         tr(iv,q) + (1-tr(iv,q-1))*tr(iv,q) + 1 - tr(iv,q) ) * tqn[iv] / 2;
      tqn[iv]  *= tr(iv,q-1);
    } 
  }

  // The tangent or ground point
  for ( iv=0; iv<nf; iv++ )    
    k(iv,0)  = ( (1-tr(iv,0))*(1+t1q[iv]*tr(iv,0)) ) * tqn[iv] / 2;
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
                    Matrix&   k,
              const Index&   start_index,
              const Index&   stop_index,
              const Matrix&   tr,
              const Index&    ground,
              const Vector&   e_ground )
{
  const Index   nf = tr.nrows(); // number of frequencies
        Index   iv;             // frequency index
        Index   q;              // LOS index (same notation as in AUG)
        Matrix   k2;             // matrix for calling other LOS WFs functions
        Vector   tr0;            // see below

  // Resize K
  k.resize( nf, start_index+1 );

  // Calculate TR0, the transmission from the platform altitude down to the
  // tangent point or the ground, and up to the platform again.
  bl( tr0, stop_index, stop_index, tr, ground, e_ground );

  // The indeces below STOP_INDEX are handled by the limb sounding function.
  sourceloswfs_limb( k2, stop_index, tr, ground, e_ground );  
  for ( iv=0; iv<nf; iv++ )
  {
    for ( q=0; q<stop_index; q++ )
      k(iv,q) = k2(iv,q);
  }

  // The indecies above STOP_INDEX are handled by the 1pass function.
  // The transmission below STOP_INDEX must here be considered.
  sourceloswfs_1pass( k2, start_index, stop_index, tr, ground, e_ground );
  for ( iv=0; iv<nf; iv++ )
  {
    for ( q=stop_index+1; q<=start_index; q++ )
      k(iv,q) = k2(iv,q)*tr0[iv];
  }

  // The platform altitude must be treated seperately
  //
  // TR0 is moved one step downwards
  q = stop_index; 
  for ( iv=0; iv<nf; iv++ )
  {
    tr0[iv] /= tr(iv,q-1)*tr(iv,q-1);    
    k(iv,q)  = ( (1-tr(iv,q-1)*tr(iv,q))*tr0[iv]*tr0[iv]*tr(iv,q-1) + 1 - 
                 tr(iv,q-1) ) / 2;
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
                    ArrayOfMatrix&   sourceloswfs,
              const Los&             los,   
              const ArrayOfMatrix&   trans,
              const Vector&          f_mono,
              const Vector&          e_ground )
{
  const Index  nza = los.start.nelem();   // number of zenith angles  

  // Resize the LOS WFs array
  sourceloswfs.resize(nza);

  // Loop zenith angles
  out3 << "    Zenith angle nr:      ";
  for (Index i=0; i<nza; i++ ) 
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    // Do something only if LOS has any points
    if ( los.p[i].nelem() > 0 )
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
                    Matrix&          k,
                    ArrayOfString&   k_names,
                    Matrix&          k_aux,
              const Los&             los,           
              const ArrayOfMatrix&   absloswfs,
              const Vector&          p_abs,
              const Vector&          t_abs,             
              const TagGroups&       tags,
              const ArrayOfMatrix&   abs_per_tg,
              const Matrix&          vmrs,
              const Vector&          k_grid,
              const ArrayOfIndex&    tg_nr,
              const String&          unit )
{
  check_lengths( p_abs, "p_abs", t_abs, "t_abs" );  
  check_length_ncol( p_abs, "p_abs", vmrs, "vmrs" );
  if ( tags.nelem() != abs_per_tg.nelem() )
    throw runtime_error(
                       "Lengths of *wfs_tgs* and *abs_per_tg* do not match." );
  if ( los.p.nelem() != absloswfs.nelem() )
    throw runtime_error(
     "The number of zenith angles is not the same in *los* and *absloswfs*." );

  // Main sizes
  const Index  nza = los.start.nelem();      // number of zenith angles  
  const Index  nv  = abs_per_tg[0].nrows();  // number of frequencies
  const Index  ntg = tg_nr.nelem();          // number of retrieval tags to do
  const Index  np  = k_grid.nelem();         // number of retrieval altitudes

  // -log(p) is used as altitude variable. The minus is included to get
  // increasing values, a demand for the grid functions. 
  Vector  lgrid;                      // -log of the retrieval pressures
  p2grid( lgrid, k_grid );
  Vector  lplos;                      // -log of the LOS pressures

  // Indices
  // IP0 and IF0 are the index off-sets for the total K matrix
        Index  itg;                       // Tag index
        Index  iza;                       // Zenith angle index
        Index  ip, ip0=0;                 // Retrieval point indices
        Index  iv, iv0;                   // Frequency indices
        Index  i1, iw;                    // weight indices

  // Other variables
        Matrix  abs;                       // absorption at k_grid
        Matrix  is;                        // matrix for storing LOS index
        Vector  w;                         // weights for LOS WFs
        Vector  a(nv);                     // temporary vector
        Index  tg;                        // present tag nr
        Vector  vmr, p, t;                 // for conversion to VMR and ND 
        Vector  nd;                        // number density

  // Set up K and additional data. Set all values of K to 0
  k.resize(nza*nv,ntg*np);
  k = 0.0;			// Matpack can set all elements like this.
  k_names.resize(ntg);
  k_aux.resize(ntg*np,2);

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
    abs.resize( nv, np );
    interpp( abs, p_abs, abs_per_tg[tg], k_grid );

    // Fill K_NAMES
    {
      ostringstream os;
      os << "Species: " << tags[tg][0].Name();
      k_names[itg] = os.str();
    }

    // Fill column 0 of K_AUX
    for ( ip=0; ip<np; ip++ )
       k_aux(ip0+ip,0) = k_grid[ip];

    // Fill column 1 of K_AUX
    //   frac : fractions of the linearisation profile
    //   vmr  : VMR
    //   nd   : number density
    //
    vmr.resize( np );
    interpp( vmr, p_abs, vmrs(tg,Range(joker)), k_grid ); 
    // The Range(joker) expression selects the entire row.

    if ( unit == "frac" )
    {
      for ( ip=0; ip<np; ip++ )
        k_aux(ip0+ip,1) = 1.0;
    }
    else if ( unit == "vmr" )
    {
      for ( ip=0; ip<np; ip++ )
      {
        for ( iv=0; iv<nv; iv++ )
          abs(iv,ip)   /= vmr[ip];
        k_aux(ip0+ip,1) = vmr[ip];
      }
    }  
    else if ( unit == "nd" )
    {
      nd.resize( np );
      interpp(  nd, p_abs, number_density(p_abs,t_abs), k_grid );
      nd *= vmr;		// Matpack can do element-vise
				// multiplication like this.
      for ( ip=0; ip<np; ip++ )
      {
        for ( iv=0; iv<nv; iv++ )
          abs(iv,ip)   /= nd[ip];
        k_aux(ip0+ip,1) = nd[ip];
      }
    }
    else
      throw runtime_error(
        "Allowed species retrieval units are \"frac\", \"vmr\" and \"nd\"."); 

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
      if ( los.p[iza].nelem() > 0 )
      {
        // Get the LOS points affected by each retrieval point
	//        lplos = -1.0 * log(los.p[iza]);
        p2grid( lplos, los.p[iza] );
        grid2grid_index( is, lplos, lgrid );

        // Loop retrieval points
        for ( ip=0; ip<np; ip++ ) 
        {
          // Check if there is anything to do
          if ( is(ip,0) >= 0 )
          {
            // Get the weights for the LOS points
	    grid2grid_weights( w, lplos, (Index)is(ip,0), (Index)is(ip,1), 
                                                                 lgrid, ip );

            // Calculate the WFs.
            // A is di/dkappa*dkappa/dkp in a compact form.
            // This is possible as the columns of dkappa/dkp are identical.  

            a = 0.0;		// Matpack can set all elements like this.

            i1 = (Index)is(ip,0);       // first LOS point to consider
            for ( iv=0; iv<nv; iv++ )
	    {
              for ( iw=i1; iw<=(Index)is(ip,1); iw++ )
                a[iv] += absloswfs[iza](iv,iw) * w[iw-i1];
              k(iv0+iv,ip0+ip) = a[iv] * abs(iv,ip);                    
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
   points (order+1) that are evenly spread between the given frequencies.
   The weighting functions are zero outside these limits.

   \retval   k            weighting function matrix
   \retval   k_names      identity name(s)
   \retval   k_aux        additional data
   \param    los          line of sight structure
   \param    absloswfs    absorption LOS Wfs
   \param    f_mono       frequency absoprtion grid
   \param    k_grid       retrieval grid
   \param    order        polynomial order
   \param    flow         lower frequency limit of fit
   \param    fhigh        upper frequency limit of fit

   \author Patrick Eriksson
   \date   2000-09-15

   Adapted to MTL. 
   \date   2001-01-05
   \author Stefan Buehler

   Made it possible to have arbitrary frequency limits.
   \date   2001-01-21
   \author Patrick Eriksson
*/
void k_contabs (
                    Matrix&          k,
                    ArrayOfString&   k_names,
                    Matrix&          k_aux,
              const Los&             los,           
              const ArrayOfMatrix&   absloswfs,
              const Vector&          f_mono,
              const Vector&          k_grid,
              const Index&           order,
              const Numeric&         flow,
              const Numeric&         fhigh )
{
  if ( los.p.nelem() != absloswfs.nelem() )
    throw runtime_error(
     "The number of zenith angles is not the same in *los* and *absloswfs*." );
  check_length_nrow( f_mono, "f_mono", absloswfs[0], 
                                                 "the matrices of absloswfs" );
  // Main sizes
  const Index  nza = los.start.nelem();     // number of zenith angles  
  const Index  np  = k_grid.nelem();        // number of retrieval altitudes
  const Index  npoints = order+1;          // number of off-set points
  const Index  nv  = f_mono.nelem();        // number of frequencies

  // Check given frequency limits
  assert( flow >= 0 );
  assert( fhigh >= 0 );
  if ( flow >= fhigh )
    throw runtime_error(
             "The lower frequency limit equals or is above the upper limit." );
  if ( flow >= f_mono[nv-1] )
    throw runtime_error(
                  "The lower frequency limit is above all values of f_mono." );
  if ( fhigh <= f_mono[0] )
    throw runtime_error(
                  "The upper frequency limit is below all values of f_mono." );

  // -log(p) is used as altitude variable. The minus is included to get
  // increasing values, a demand for the grid functions. 
  Vector  lgrid;
  p2grid( lgrid, k_grid );
  Vector  lplos;                      // -log of the LOS pressures

  // Indices
  // IP0 and IF0 are the index off-sets for the total K matrix
        Index  ipoint;                    // Off-set point index
        Index  iza;                       // Zenith angle index
        Index  ip, ip0=0;                 // Retrieval point indices
        Index  iv, iv0;                   // Frequency indices
        Index  i1, iw;                    // weight indices

  // Determine first and last frequency index inside given limits
  Index   ilow, ihigh;
  for( ilow=0; ilow<nv && f_mono[ilow] < flow; ilow++ )
    {}
  for( ihigh=ilow; ihigh<nv && f_mono[ihigh] <= fhigh; ihigh++ )
    {}

  // Other variables
        Vector  fpoints;                   // frequencies of the off-set points
        Vector  b(nv);                     // fit base function
        Matrix  is;                        // matrix for storing LOS index
        Vector  w;                         // weights for LOS WFs
        Vector  a(nv);                     // temporary vector

  // Check that the selected polynomial order
  if ( order < 0 )
    throw runtime_error("The polynomial order must be >= 0."); 

  // Set up K and additional data. Set all values of K to 0
  k.resize(nza*nv,npoints*np);
  k = 0.0;			// Matpack can set all elements like this.
  k_names.resize(npoints);
  k_aux.resize(npoints*np,2);

  // Calculate the frequencies of the off-set points
  if ( npoints > 1 )
    nlinspace( fpoints, f_mono[ilow], f_mono[ihigh], npoints );
  else
  {
    fpoints.resize( 1 );
    fpoints[0] = ( f_mono[ilow] + f_mono[ihigh] ) / 2.0;
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
       k_aux(ip0+ip,0) = k_grid[ip];
       k_aux(ip0+ip,1) = 0.0;
    }

    // Set-up base vector for the present fit point 
    b = 1.0;			// Matpack can set all elements like this.
    if ( npoints > 1 )
    {
      for ( ip=0; ip<npoints; ip++ )
      {
        if ( ip != ipoint )
	{ 
          for ( iv=ilow; iv<=ihigh; iv++ )
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
      if ( los.p[iza].nelem() > 0 )
      {
        // Get the LOS points affected by each retrieval point
	//        lplos = -1.0 * log(los.p[iza]);
        p2grid( lplos, los.p[iza] );
        grid2grid_index( is, lplos, lgrid );

        // Loop retrieval points
        for ( ip=0; ip<np; ip++ ) 
        {
          // Check if there is anything to do
          if ( is(ip,0) >= 0 )
          {
            // Get the weights for the LOS points
	    grid2grid_weights( w, lplos, (Index)is(ip,0), (Index) is(ip,1),
                                                                 lgrid, ip );

            // Calculate the WFs.
            // A is di/dkappa*dkappa/dkp in a compact form.
            // This is possible as the columns of dkappa/dkp are identical.  

            a = 0.0;		// Matpack can set all elements like this.

            i1 = (Index)is(ip,0);       // first LOS point to consider
            for ( iv=ilow; iv<=ihigh; iv++ )
	    {
              for ( iw=i1; iw<=(Index)is(ip,1); iw++ )
                a[iv] += absloswfs[iza](iv,iw) * w[iw-i1];
              k(iv0+iv,ip0+ip) = a[iv] * b[iv];                    
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
              Matrix&                     k,
              ArrayOfString&              k_names,
              Matrix&                     k_aux,
        const TagGroups&                  tag_groups,
        const Los&                        los,           
        const ArrayOfMatrix&              absloswfs,
        const Vector&                     f_mono,
        const Vector&                     p_abs,
        const Vector&                     t_abs,
        const Vector&                     n2_abs,	   
        const Vector&                     h2o_abs,	   
        const Matrix&                     vmrs,
        const ArrayOfArrayOfLineRecord&   lines_per_tg,
        const ArrayOfLineshapeSpec&       lineshape,
        const Matrix&                     abs,            
        const ArrayOfMatrix&              trans,
        const Vector&                     e_ground,
        const Vector&                     k_grid,
        const ArrayOfString&              cont_description_names,
	const ArrayOfVector& 	          cont_description_parameters,
        const ArrayOfString&              cont_description_models )
{
  // Main sizes
  const Index  nza = los.start.nelem();     // number of zenith angles  
  const Index  nv  = f_mono.nelem();        // number of frequencies
  const Index  np  = k_grid.nelem();        // number of retrieval altitudes

  // -log(p) is used as altitude variable. The minus is included to get
  // increasing values, a demand for the grid functions. 
  //  const Vector  lgrid=-1.0*log(k_grid);
  Vector  lgrid;
  p2grid( lgrid, k_grid );
  Vector  lplos;                     // -log of the LOS pressures

  // Indices
  // IP0 and IF0 are the index off-sets for the total K matrix
        Index  iza;                       // zenith angle index
        Index  ip;                        // retrieval point index
        Index  iv, iv0=0;                 // frequency indices
        Index  i1, iw;                    // weight indices

  // Other variables
        Vector  t(k_grid.nelem());          // temperature at retrieval points
        Matrix  abs1k;                     // absorption for t_abs+1K
        Matrix  dabs_dt;                   // see below
 ArrayOfMatrix  abs_dummy;                 // dummy absorption array
 ArrayOfMatrix  sloswfs;                   // source LOS WFs
	Matrix  is;                        // matrix for storing LOS index
        Vector  w;                         // weights for LOS WFs
        Vector  a(nv), b(nv), pl(f_mono.nelem());  // temporary vectors

 // The scalars are declared to be double to avoid possible numerical problems
 // when using float
        double  c,d;                       // temporary values


  // Set up K and additional data. Set all values of K to 0
  k.resize(nza*nv,np);
  k = 0.0;			// Matpack can set all elements like this.
  k_names.resize(1);
  k_names[0] = "Temperature: no hydrostatic eq.";
  k_aux.resize(np,2);
  interpp( t, p_abs, t_abs, k_grid ); 
  for ( ip=0; ip<np; ip++ )
  {
     k_aux(ip,0) = k_grid[ip];
     k_aux(ip,1) = t[ip];
  }

  // Calculate absorption for t_abs + 1K to estimate the temperature derivative
  // dabs/dt, the temperature derivative of the absorption at k_grid
  out2 << "  Calculating absorption for t_abs + 1K\n";
  out2 << "  ----- Messages from absCalc: -----\n";
  //
  {
    // Dummy should hold t_abs + 1:
    Vector dummy(t_abs);	// Matpack can initialize a Vector
				// from another Vector. 
    dummy += 1;			// Matpack can add element-vise like this.

    absCalc( abs1k, abs_dummy, tag_groups, f_mono, p_abs, dummy, n2_abs, 
             h2o_abs, vmrs, lines_per_tg, lineshape, 
	     cont_description_names, 
	     cont_description_models,
             cont_description_parameters);
  }
  abs_dummy.resize(0);
  //
  out2 << "  ----- Back from absCalc ----------\n";
  //
  // Compute abs1k = abs1k - abs:
  abs1k -= abs;			// Matpack can subtract element-vise like this.

  dabs_dt.resize( abs1k.nrows(), k_grid.nelem() );
  interpp( dabs_dt, p_abs, abs1k, k_grid );
  abs1k.resize(0,0);

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
    if ( los.p[iza].nelem() > 0 )
    {
      // Get the LOS points affected by each retrieval point
      //      lplos = -1.0 * log(los.p[iza]);
      p2grid( lplos, los.p[iza] );
      grid2grid_index( is, lplos, lgrid );

      // Loop retrieval points
      for ( ip=0; ip<np; ip++ ) 
      {
        // Check if there is anything to do
        if ( is(ip,0) >= 0 )
        {
          // Get the weights for the LOS points
	  grid2grid_weights( w, lplos, (Index)is(ip,0), (Index) is(ip,1), 
                                                                 lgrid, ip );

          // Calculate the WFs.
          // A is di/dsigma*dsigma/dSp in a compact form.
          // B is di/dkappa*dkappa/dkp in a compact form.
          // This is possible as the columns of dkappa/dkp are identical and  
	  // that dkappa/dkp = dsigma/dSp
          // C is just a temporary value
	  // PL is the Planck function for the present temperature value
	  //
          a = 0.0;		// Matpack can set all elements like this.
	  b = 0.0;                    
          c  = PLANCK_CONST / BOLTZMAN_CONST / t[ip];
          planck( pl, f_mono, t[ip] );
          i1 = (Index)is(ip,0);       // first LOS point to consider
          //
          for ( iv=0; iv<nv; iv++ )
	  {
            for ( iw=i1; iw<=(Index)is(ip,1); iw++ )
	    {
              a[iv] += sloswfs[iza](iv,iw) * w[iw-i1];
              b[iv] += absloswfs[iza](iv,iw) * w[iw-i1];
	    }
            d = c * f_mono[iv];
            k(iv0+iv,ip) = a[iv] * d/t[ip] / (1-exp(-d)) * pl[iv] +
                                                      b[iv] * dabs_dt(iv,ip);
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
		   const ArrayOfString& tags)
{
  wfs_tag_groups = TagGroups(tags.nelem());

  // Each element of the array of Strings tags defines one tag
  // group. Let's work through them one by one.
  for ( Index i=0; i<tags.nelem(); ++i )
    {
      // There can be a comma separated list of tag definitions, so we
      // need to break the String apart at the commas.
      ArrayOfString tag_def;

      bool go_on = true;
      String these_tags = tags[i];
      while (go_on)
	{
	  Index n = these_tags.find(',');
	  if ( n == these_tags.npos ) // npos indicates `not found'
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
      // the outer Array of an Array<Array>> pointing to the same data
      // after creation. So we need to fix this explicitly:
      wfs_tag_groups[i] = Array<OneTag>();

      for ( Index s=0; s<tag_def.nelem(); ++s )
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
  for ( Index i=0; i<wfs_tag_groups.nelem(); ++i )
    {
      out3 << "\n  " << i << ":";
      for ( Index s=0; s<wfs_tag_groups[i].nelem(); ++s )
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
                    ArrayOfMatrix&   absloswfs,
	      const Index&             emission,
              const Los&             los,   
              const ArrayOfMatrix&   source,
              const ArrayOfMatrix&   trans,
              const Vector&          y,
              const Vector&          y_space,
              const Vector&          f_mono,
              const Vector&          e_ground,
              const Numeric&         t_ground )
{
  check_if_bool( emission, "emission" );                                      

  if ( emission == 0 )
    absloswfs_tau ( absloswfs, los, f_mono );
  else
    absloswfs_rte ( absloswfs, los, source, trans, y, y_space, f_mono,     
                                                          e_ground, t_ground );

  // The equations used can produce NaNs and Infs when the optical thickness
  // is very high. This corresponds to extremly low values for the ABS LOS WFs
  // so we set NaNs and Infs to zero.
  Index     irow, icol, ncol;
  Numeric   w;
  for( Index im=0; im<absloswfs.nelem(); im++ )
  {
    for( irow=0; irow<absloswfs[im].nrows(); irow++ )
    {
      ncol = absloswfs[im].ncols();
      for( icol=0; icol<ncol; icol++ )
      {
        w = absloswfs[im](irow,icol);
        if ( isnan(w) || isinf(w) )
          absloswfs[im](irow,icol) = 0.0;
      }
    }
  }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kSpecies (
                    Matrix&          k,
                    ArrayOfString&   k_names,
                    Matrix&          k_aux,
              const Los&             los,           
              const ArrayOfMatrix&   absloswfs,
              const Vector&          p_abs,
              const Vector&          t_abs,             
              const TagGroups&       wfs_tgs,
              const ArrayOfMatrix&   abs_per_tg,
              const Matrix&          vmrs,
              const Vector&          k_grid,
              const String&          unit )
{
  // Check of input is performed in k_species

  const Index   ntg = wfs_tgs.nelem();     // number of retrieval tags
  ArrayOfIndex  tg_nr(ntg);

  for ( Index i=0; i<ntg; i++ )
    tg_nr[i] = i;

  k_species( k, k_names, k_aux, los, absloswfs, p_abs, t_abs, 
	     wfs_tgs, abs_per_tg, vmrs, k_grid, tg_nr, unit );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kSpeciesSingle (
                    Matrix&          k,
                    ArrayOfString&   k_names,
                    Matrix&          k_aux,
              const Los&             los,           
              const ArrayOfMatrix&   absloswfs,
              const Vector&          p_abs,
              const Vector&          t_abs,             
              const TagGroups&       wfs_tgs,
              const ArrayOfMatrix&   abs_per_tg,
              const Matrix&          vmrs,
              const Vector&          k_grid,
              const String&          tg,
              const String&          unit )
{
  // Check of input is performed in k_species

  ArrayOfString  tg_name(1);
  tg_name[0] = tg;

  ArrayOfIndex   tg_nr; 
  get_tagindex_for_Strings( tg_nr, wfs_tgs, tg_name );
  
  k_species( k, k_names, k_aux, los, absloswfs, p_abs, t_abs, 
                 	     wfs_tgs, abs_per_tg, vmrs, k_grid, tg_nr, unit );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-01-21
*/
void kContAbs (
                    Matrix&          k,
                    ArrayOfString&   k_names,
                    Matrix&          k_aux,
              const Los&             los,           
              const ArrayOfMatrix&   absloswfs,
              const Vector&          f_mono,
              const Vector&          k_grid,
	      const Index&           order,
              const Numeric&         f_low,
              const Numeric&         f_high )
{
  // Input is checked in k_contabs  

  Numeric f1=f_low, f2=f_high;

  if ( f1 < 0 )
    f1 = first( f_mono );
  if ( f2 < 0 )
    f2 = last( f_mono );

  k_contabs( k, k_names, k_aux, los, absloswfs, f_mono, k_grid, order, f1, f2);
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-04-18
*/
void kTemp (
                Matrix&                      k,
                ArrayOfString&               k_names,
                Matrix&                      k_aux,
          const TagGroups&                   tgs,
          const Vector&                      f_mono,
          const Vector&                      p_abs,
          const Vector&                      t_abs,
          const Vector&                      n2_abs,
          const Vector&                      h2o_abs,
          const Matrix&                      vmrs,
	  const Matrix&                      abs0,
          const ArrayOfArrayOfLineRecord&    lines_per_tg,
          const ArrayOfLineshapeSpec&        lineshape,
          const Vector&                      e_ground,
          const Index&                       emission,
          const Vector&                      k_grid,
          const ArrayOfString&               cont_description_names,
          const ArrayOfVector& 	             cont_description_parameters,
          const ArrayOfString&               cont_description_models,
          const Los&                         los,           
          const ArrayOfMatrix&               absloswfs,
          const ArrayOfMatrix&               trans,
    	  const Numeric&    		     z_plat,
    	  const Vector&     		     za,
    	  const Numeric&    		     l_step,
    	  const Vector&     		     z_abs,
    	  const Index&        		     refr,
    	  const Index&        		     refr_lfac,
    	  const Vector&     		     refr_index,
    	  const Numeric&    		     z_ground,
          const Numeric&                     t_ground,
    	  const Vector&     		     y_space,
    	  const Numeric&    		     r_geoid,
          const Vector&     		     hse,
	  // Keywords
          const Index&                       kw_hse,
          const Index&                       kw_fast )
{
  check_if_bool( kw_hse, "hse keyword" );
  check_if_bool( kw_fast, "fast keyword" );      
  check_if_bool( emission, "emission" );
  check_if_bool( refr, "refr" );
  check_lengths( p_abs, "p_abs", t_abs, "t_abs" );  
  check_lengths( p_abs, "p_abs", z_abs, "z_abs" );  
  check_lengths( p_abs, "p_abs", h2o_abs, "h2o_abs" );  
  check_lengths( p_abs, "p_abs", n2_abs, "n2_abs" );  
  check_length_ncol( p_abs, "p_abs", abs0, "abs" );
  check_length_ncol( p_abs, "p_abs", vmrs, "vmrs" );
  check_length_nrow( f_mono, "f_mono", abs0, "abs" );
  if ( los.p.nelem() != za.nelem() )
    throw runtime_error(
               "The number of zenith angles of *za* and *los* are different.");
  //
  if( refr ) 
    check_lengths( p_abs, "p_abs", refr_index, "refr_index" );  
  //
  if ( !kw_hse )
  {
    if ( los.p.nelem() != trans.nelem() )
      throw runtime_error(
            "The number of zenith angles of *los* and *trans* are different.");
    check_length_nrow( f_mono, "f_mono", trans[0], 
                                       "the transmission matrices (in trans)");
    if ( los.p.nelem() != absloswfs.nelem() )
      throw runtime_error(
      "The number of zenith angles is not the same in *los* and *absloswfs*.");
  }
  else
  { 
    if ( hse[0]==0 )
      throw runtime_error( "Hydrostatic eq. must be considered generally" 
                           "when calculating WFs with hydrostatic eq.");
  }
  //
  if ( any_ground(los.ground) )  
  {
    if ( t_ground <= 0 )
      throw runtime_error(
          "There are intersections with the ground, but the ground\n"
          "temperature is set to be <=0 (are dummy values used?).");
    if ( e_ground.nelem() != f_mono.nelem() )
      throw runtime_error(
          "There are intersections with the ground, but the frequency and\n"
          "ground emission vectors have different lengths (are dummy values\n"
          "used?).");
  }
  if ( emission ) 
    check_lengths( f_mono, "f_mono", y_space, "y_space" );

  //
  // Three options:
  // 1. no hydrostatic eq., use analytical expressions
  // 2. hydrostatic eq., fast version
  // 3. hydrostatic eq., accurate version
  //


  // No hydrostatic eq., use analytical expressions
  //---------------------------------------------------------------------------
  if ( !kw_hse )
  {
    if ( !emission )
      throw runtime_error(
          "Analytical expressions for temperature and no emission have not\n"
          "yet been implemented. Sorry!");

    k_temp_nohydro( k, k_names, k_aux, tgs, los, absloswfs, f_mono, 
     p_abs, t_abs, n2_abs, h2o_abs, vmrs, lines_per_tg, lineshape, abs0, trans,
     e_ground, k_grid, cont_description_names, cont_description_parameters,
     cont_description_models);
  }


  // Hydrostatic eq., fast version
  //---------------------------------------------------------------------------
  else if ( kw_fast )
  {
    // Main sizes
    const Index  nza  = za.nelem();          // number of zenith angles  
    const Index  nv   = f_mono.nelem();      // number of frequencies
    const Index  np   = k_grid.nelem();      // number of retrieval altitudes
    const Index  nabs = p_abs.nelem();       // number of absorption altitudes
  
    // Vectors for the reference state
    Vector z0(nabs), y0, t0(np);
  
    // Local copy of hse
    Vector hse_local( hse );
  
    // Calculate reference z_abs with a high number of iterations
    hse_local[4] = 5;
    z0 = z_abs;
    hseCalc( z0, p_abs, t_abs, h2o_abs, r_geoid,  hse_local );
    hse_local[4] = hse[4];
  
    // Calculate absorption for + 1K
    Matrix         abs1k, abs(nv,nabs);
    ArrayOfMatrix  abs_dummy;
    //
    {
      Vector  t(t_abs);
      t += 1.;
      //
      out1 << "  Calculating absorption for t_abs + 1K \n";
      out2 << "  ----- Messages from absCalc: --------\n";
      absCalc( abs1k, abs_dummy, tgs, f_mono, p_abs, t, n2_abs, h2o_abs, vmrs, 
                 lines_per_tg, lineshape, 
	       cont_description_names, 
	       cont_description_models,
               cont_description_parameters);
    }
    // Calculate reference spectrum
    out1 << "  Calculating reference spectrum\n";
    out2 << "  ----- Messages from losCalc: --------\n";
    Los    los;
    Vector z_tan;
    losCalc( los, z_tan, z_plat, za, l_step, p_abs, z_abs, refr, refr_lfac, 
					       refr_index, z_ground, r_geoid );
    out2 << "  -------------------------------------\n";
    out2 << "  ----- Messages from sourceCalc: -----\n";
    ArrayOfMatrix source, trans;
    sourceCalc( source, emission, los, p_abs, t_abs, f_mono );
    out2 << "  -------------------------------------\n";
    out2 << "  ----- Messages from transCalc: ------\n";
    transCalc( trans, los, p_abs, abs0 );
    out2 << "  -------------------------------------\n";
    out2 << "  ----- Messages from yRte: -----------\n";
    yCalc( y0, emission, los, f_mono, y_space, source, trans, 
							  e_ground, t_ground );
    out2 << "  -------------------------------------\n";
  
    // Allocate K and fill aux. variables
    k.resize(nza*nv,np);
    k_names.resize(1);
    k_names[0] = "Temperature: with hydrostatic eq.";
    k_aux.resize(np,2);
    interpp( t0, p_abs, t_abs, k_grid ); 
    for ( Index ip=0; ip<np; ip++ )
    {
       k_aux(ip,0) = k_grid[ip];
       k_aux(ip,1) = t0[ip];
    }
  
    // Determine conversion between grids        
    Matrix is;
    Vector lpabs, lgrid;
    p2grid( lpabs, p_abs );
    p2grid( lgrid, k_grid );
    grid2grid_index( is, lpabs, lgrid );
  
    // Loop retrieval altitudes and calculate new spectra
    //
    Vector y, w;
    Index  i1, iw, iv;
    //
    for ( Index ip=0; ip<np; ip++ )
    {
      out1 << "  Doing altitude " << ip+1 << "/" << np << "\n";   
  
      // Create absorption matrix corresponding to temperature disturbance
      grid2grid_weights( w, lpabs, Index(is(ip,0)), Index(is(ip,1)), 
								   lgrid, ip );
      i1 = Index( is(ip,0) );    // first p_abs point to consider
      abs = abs0;	

      for ( iw=i1; iw<=Index(is(ip,1)); iw++ )
      {
	for ( iv=0; iv<nv; iv++ )
	  abs(iv,iw) = (1-w[iw-i1])*abs0(iv,iw) + w[iw-i1]*abs1k(iv,iw);
      }
  
      out2 << "  ----- Messages from losCalc: --------\n";
      losCalc( los, z_tan, z_plat, za, l_step, p_abs, z_abs, refr, refr_lfac, 
					       refr_index, z_ground, r_geoid );
      out2 << "  -------------------------------------\n";
      out2 << "  ----- Messages from sourceCalc: -----\n";
      ArrayOfMatrix source, trans;
      sourceCalc( source, emission, los, p_abs, t_abs, f_mono );
      out2 << "  -------------------------------------\n";
      out2 << "  ----- Messages from transCalc: ------\n";
      transCalc( trans, los, p_abs, abs );
      out2 << "  -------------------------------------\n";
      out2 << "  ----- Messages from yRte: -----------\n";
      yCalc( y, emission, los, f_mono, y_space, source, trans, 
							  e_ground, t_ground );
      out2 << "  -------------------------------------\n";
  
      // Fill K
      for ( iv=0; iv<nza*nv; iv++ )
	k(iv,ip) = y[iv] - y0[iv];
    }
  }


  // Hydrostatic eq., accurate version
  //---------------------------------------------------------------------------
  else
  {
    // Main sizes
    const Index  nza  = za.nelem();          // number of zenith angles  
    const Index  nv   = f_mono.nelem();      // number of frequencies
    const Index  np   = k_grid.nelem();      // number of retrieval altitudes
    const Index  nabs = p_abs.nelem();       // number of absorption altitudes
  
    // Vectors for the reference state
    Vector z0(nabs), y0, t0(np);
  
    // Local copy of hse
    Vector hse_local( hse );
  
    // Calculate reference z_abs with a high number of iterations
    hse_local[4] = 5;
    z0 = z_abs;	
    hseCalc( z0, p_abs, t_abs, h2o_abs, r_geoid,  hse_local );
    hse_local[4] = hse[4];
  
    // Calculate reference spectrum
    out1 << "  Calculating reference spectrum\n";
    out2 << "  ----- Messages from losCalc: --------\n";
    Los    los;
    Vector z_tan;
    losCalc( los, z_tan, z_plat, za, l_step, p_abs, z_abs, refr, refr_lfac, 
	        			       refr_index, z_ground, r_geoid );
    out2 << "  -------------------------------------\n";
    out2 << "  ----- Messages from sourceCalc: -----\n";
    ArrayOfMatrix source, trans;
    sourceCalc( source, emission, los, p_abs, t_abs, f_mono );
    out2 << "  -------------------------------------\n";
    out2 << "  ----- Messages from transCalc: ------\n";
    transCalc( trans, los, p_abs, abs0 );
    out2 << "  -------------------------------------\n";
    out2 << "  ----- Messages from yRte: -----------\n";
    yCalc( y0, emission, los, f_mono, y_space, source, trans, 
							  e_ground, t_ground );
    out2 << "  -------------------------------------\n";
  
    // Allocate K and fill aux. variables
    k.resize(nza*nv,np);
    k_names.resize(1);
    k_names[0] = "Temperature: with hydrostatic eq.";
    k_aux.resize(np,2);
    interpp( t0, p_abs, t_abs, k_grid ); 
    for ( Index ip=0; ip<np; ip++ )
    {
       k_aux(ip,0) = k_grid[ip];
       k_aux(ip,1) = t0[ip];
    }
  
    // Determine conversion between grids        
    Matrix is;
    Vector lpabs, lgrid;
    p2grid( lpabs, p_abs );
    p2grid( lgrid, k_grid );
    grid2grid_index( is, lpabs, lgrid );
  
    // Loop retrieval altitudes and calculate new spectra
    //
    Matrix         abs;
    ArrayOfMatrix  abs_dummy;
    Vector y, t(nabs), w;
    Index  i1, iw, iv;
    //
    for ( Index ip=0; ip<np; ip++ )
    {
      out1 << "  Doing altitude " << ip+1 << "/" << np << "\n";   
  
      // Create disturbed temperature profile
      grid2grid_weights( w, lpabs, Index(is(ip,0)), Index(is(ip,1)), 
								   lgrid, ip );
      i1 = Index( is(ip,0) );    // first p_abs point to consider
      t = t_abs;			// Matpack can copy the contents of
					// vectors like this. The dimensions
					// must be the same! 
      for ( iw=i1; iw<=Index(is(ip,1)); iw++ )
	t[iw] += w[iw-i1];
  
      out2 << "  ----- Messages from absCalc: --------\n";
      absCalc( abs, abs_dummy, tgs, f_mono, p_abs, t, n2_abs, h2o_abs, vmrs, 
	       lines_per_tg, lineshape, 
	       cont_description_names, 
	       cont_description_models,
               cont_description_parameters);
      out2 << "  ----- Messages from losCalc: --------\n";
      losCalc( los, z_tan, z_plat, za, l_step, p_abs, z_abs, refr, refr_lfac, 
					       refr_index, z_ground, r_geoid );
      out2 << "  -------------------------------------\n";
      out2 << "  ----- Messages from sourceCalc: -----\n";
      ArrayOfMatrix source, trans;
      sourceCalc( source, emission, los, p_abs, t_abs, f_mono );
      out2 << "  -------------------------------------\n";
      out2 << "  ----- Messages from transCalc: ------\n";
      transCalc( trans, los, p_abs, abs );
      out2 << "  -------------------------------------\n";
      out2 << "  ----- Messages from yRte: -----------\n";
      yCalc( y, emission, los, f_mono, y_space, source, trans, 
							  e_ground, t_ground );
      out2 << "  -------------------------------------\n";
  
      // Fill K
      for ( iv=0; iv<nza*nv; iv++ )
	k(iv,ip) = y[iv] - y0[iv];
    }
  }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-21
 */
void kPointingOffSet(
                    Matrix&          k,
                    ArrayOfString&   k_names,
                    Matrix&          k_aux,
              const Numeric&         z_plat,
              const Vector&          za_pencil,
              const Numeric&         l_step,
              const Vector&          p_abs,
              const Vector&          z_abs,
              const Vector&          t_abs,
              const Vector&          f_mono,
              const Index&           refr,
              const Index&           refr_lfac,
              const Vector&          refr_index,
              const Numeric&         z_ground,
              const Numeric&         r_geoid,
              const Matrix&          abs,
  	      const Index&           emission,
              const Vector&          y_space,
              const Vector&          e_ground,
              const Numeric&         t_ground,
              const Vector&          y,
              const Numeric&         delta )
{
  check_if_bool( emission, "emission" );
  check_if_bool( refr, "refr" );
  check_lengths( p_abs, "p_abs", t_abs, "t_abs" );  
  check_lengths( p_abs, "p_abs", z_abs, "z_abs" );  
  check_length_ncol( p_abs, "p_abs", abs, "abs" );
  check_length_nrow( f_mono, "f_mono", abs, "abs" );
  if ( emission ) 
    check_lengths( f_mono, "f_mono", y_space, "y_space" );
  if ( refr )
    check_lengths( p_abs, "p_abs", refr_index, "refr_index" );  


  // Create new zenith angle grid
  //  const Index  nza = za_pencil.nelem();
  Vector za_new( za_pencil );	// Matpack can initialize a
				// new Vector from another
				// Vector.
  za_new += delta;		// With Matpack you can add delta to all
				// Vector elements like this.

  out2 << "  ----- Messages from losCalc: --------\n";
  Los    los;
  Vector z_tan;
  losCalc( los, z_tan, z_plat, za_new, l_step, p_abs, z_abs, refr, refr_lfac, 
                                              refr_index, z_ground, r_geoid );
  out2 << "  -------------------------------------\n";

  ArrayOfMatrix source, trans;
  out2 << "  ----- Messages from sourceCalc: -----\n";
  sourceCalc( source, emission, los, p_abs, t_abs, f_mono );
  out2 << "  -------------------------------------\n";
  out2 << "  ----- Messages from transCalc: ------\n";
  transCalc( trans, los, p_abs, abs );
  out2 << "  -------------------------------------\n";
  Vector y_new;
  out2 << "  ----- Messages from yRte: -----------\n";
  yCalc( y_new, emission, los, f_mono, y_space, source, trans, 
                                                         e_ground, t_ground );
  out2 << "  -------------------------------------\n";

  // Make k one-column matrix of the right size:
  const Index   nv = y.nelem();
  k.resize( nv, 1 );

  // k = (y_new - y) / delta
  k = y_new;			
  k -= y;
  k /= delta;

  /* With Matpack, you can use scalar, vector, or matrix +=,-=,*=, and
     /= operators. They all act element-vise. 
     Vectors behave like 1-column matrices, therefore copying y_new to
     k works. 

     This is the old code. Isn't it much nicer now?
     Vector dummy( nv );
     add( y_new, scaled(y,-1), dummy );
     copy( scaled(dummy,1/delta), columns(k)[0] ); */

  k_names.resize( 1 );
  k_names[0] = "Pointing: off-set";
  k_aux.resize( 1, 2 );
  k_aux = 0.0;			// Matpack can set all elements like this.
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-09-18
 */
void kEground(
                    Matrix&          k,
                    ArrayOfString&   k_names,
                    Matrix&          k_aux,
              const Vector&          za_pencil,
              const Vector&          f_mono,
  	      const Index&           emission,
              const Vector&          y_space,
              const Vector&          e_ground,
              const Numeric&         t_ground,
              const Los&             los,           
              const ArrayOfMatrix&   source,
              const ArrayOfMatrix&   trans,
              const Index&           single_e )
{
  check_if_bool( emission, "emission" );                                      
  check_if_bool( emission, "single_e" );                                      
  check_lengths( f_mono, "f_mono", e_ground, "e_ground" );
  if ( los.p.nelem() != za_pencil.nelem() )
    throw runtime_error(
               "The number of zenith angles of *za* and *los* are different.");
  //
  if ( los.p.nelem() != trans.nelem() )
    throw runtime_error(
            "The number of zenith angles of *los* and *trans* are different.");
  check_length_nrow( f_mono, "f_mono", trans[0], 
                                       "the transmission matrices (in trans)");
  //
  if ( emission )
  {
    if ( los.p.nelem() != source.nelem() )
      throw runtime_error(
           "The number of zenith angles of *los* and *source* are different.");
    check_length_nrow( f_mono, "f_mono", source[0], 
  				            "the source matrices (in source)");
  }
  //
  if ( any_ground(los.ground) )  
  {
    if ( t_ground <= 0 )
      throw runtime_error(
          "There are intersections with the ground, but the ground\n"
          "temperature is set to be <=0 (are dummy values used?).");
    if ( e_ground.nelem() != f_mono.nelem() )
      throw runtime_error(
          "There are intersections with the ground, but the frequency and\n"
          "ground emission vectors have different lengths (are dummy values\n"
          "used?).");
  }
  //
  if ( single_e ) 
  {
    for ( INDEX iv=1; iv<f_mono.nelem(); iv++ )
    {
      if ( e_ground[iv] != e_ground[0] )
        throw runtime_error(
          "A single ground emission value is assumed, but all values of\n"
          "*e_ground* are not the same.");
    }    
  }


  // Main sizes
  const Index  nza  = za_pencil.nelem();   // number of zenith angles  
  const Index  nv   = f_mono.nelem();      // number of frequencies

  // Loop variables
  Index iv, iza;

  // Make calculations assuming a single e
  Vector ksingle( nza*nv );

  // Emission
  if ( emission )
  {
    // Set all values to zero (k is zero if no ground intersection)
    ksingle = 0.0;

    // Local vectors
    Vector bground(nv), y_in(nv), tr(nv);

    for ( iza=0; iza<nza; iza++ )
    {
      if ( los.ground[iza] )
      {
        // Calculate blackbody radiation of the ground
        planck( bground, f_mono, t_ground );

        // Calculate the intensity reflected by the ground
        rte( y_in, los.start[iza], los.ground[iza]-1, trans[iza], source[iza], 
                                               y_space, 0, e_ground, bground );

        // Calculate transmission from the ground to the sensor
        tr = 1.0;
        bl_iterate( tr, los.ground[iza]-1, los.stop[iza]-1, trans[iza], nv );

        for ( iv=0; iv<nv; iv++ )
          ksingle[iza*nv+iv] = ( bground[iv] - y_in[iv] ) * tr[iv];      
      }
    }
  }

  // Optical thicknesses
  else
  {
    Numeric a;
    for ( iv=0; iv<nv; iv++ )
    {
      a = 1 / ( 1 - e_ground[iv] );
      for ( iza=0; iza<nza; iza++ )
        ksingle[iza*nv+iv] = a;
    }
  }

  // A single value for the ground emission, k is a column vector
  if ( single_e )
  {
    k_names.resize( 1 );
    k_names[0] = "Ground emission (single)";
    //
    k.resize( nza*nv, 1 );
    k = ksingle;
    //
    k_aux.resize( 1, 2 );
    k_aux(0,0) = 0;
    k_aux(0,1) = e_ground[0];    
  }

  // An e for each monochromatic frequency, k has nv columns
  else
  {
    k_names.resize( 1 );
    k_names[0] = "Ground emission";
    //
    k.resize( nza*nv, nv );
    k = 0.0;
    k_aux.resize( nv, 2 );
    //
    for ( iv=0; iv<nv; iv++ )
    {
      for ( iza=0; iza<nza; iza++ )
        k(iza*nv+iv,iv) = ksingle[iza*nv+iv];
      k_aux(iv,0) = 0;
      k_aux(iv,1) = e_ground[iv];    
    }
  }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-21
 */
void kCalibration(
                    Matrix&          k,
                    ArrayOfString&   k_names,
                    Matrix&          k_aux,
              const Vector&          y,
              const Vector&          f_mono,
              const Vector&          y0,
       	      const String&          name )
{
  check_lengths( f_mono, "f_mono", y0, "y0" );

  const Index   ny = y.nelem();
  const Index   nf = f_mono.nelem();
  const Index   nza = Index( floor(ny/nf) );

  // Make k one-column matrix of the right size:
  k.resize( ny, 1 );

  // k = y - y0
  Index j,i,i0;
  for ( j=0; j<nza; j++ )    
  {
    i0 = j*nf;
    for ( i=0; i<nf; i++ )
      k(i0+i,0) = y[i0+i] - y0[i];
  }

  k_names.resize( 1 );
  k_names[0] = "Calibration: scaling";
  k_aux.resize( 1, 2 );
  k_aux = 0.0;		
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
                    Matrix&          k,
                    ArrayOfString&   k_names,
                    Matrix&          k_aux,
              const Vector&          y0,      
              const Vector&          y,
              const String&          name,
              const Numeric&         delta,
              const Numeric&         grid,
              const Numeric&         apriori )
{
  check_lengths( y, "y", y0, "y0" );

  // Make k one-column matrix of the right size:
  k.resize( y.nelem(), 1 );

  // Compute (y-y0)
  k = y;
  k -= y0;
  k /= delta;

  k_names.resize(1);
  k_names[0] = name;
  k_aux.resize(1,2);
  k_aux(0,0) = grid;
  k_aux(0,1) = apriori;
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kxInit (
                    Matrix&          kx,
                    ArrayOfString&   kx_names,
                    ArrayOfIndex&    kx_lengths,
                    Matrix&          kx_aux )
{
  kx.resize(         0, 0 );
  kx_names.resize(   0    );
  kx_lengths.resize( 0    );
  kx_aux.resize(     0, 0 );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kbInit (
                    Matrix&          kb,
                    ArrayOfString&   kb_names,
                    ArrayOfIndex&    kb_lengths,
                    Matrix&          kb_aux )
{
  kxInit( kb, kb_names, kb_lengths, kb_aux );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kxAppend (
		    Matrix&          kx,
		    ArrayOfString&   kx_names,
		    ArrayOfIndex&    kx_lengths,
		    Matrix&          kx_aux,
              const Matrix&          k,
              const ArrayOfString&   k_names,
              const Matrix&          k_aux )
{
  k_append( kx, kx_names, kx_lengths, kx_aux, k, k_names, k_aux );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void kbAppend (
		    Matrix&          kb,
		    ArrayOfString&   kb_names,
		    ArrayOfIndex&    kb_lengths,
		    Matrix&          kb_aux,
              const Matrix&          k,
              const ArrayOfString&   k_names,
              const Matrix&          k_aux )
{
  k_append( kb, kb_names, kb_lengths, kb_aux, k, k_names, k_aux );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-21
*/
void kxAllocate (
                    Matrix&          kx,
                    ArrayOfString&   kx_names,
                    ArrayOfIndex&    kx_lengths,
                    Matrix&          kx_aux,
              const Vector&          y,
              const String&          y_name,
              const Index&             ni,
              const Index&             nx )
{
  const Index  ny = y.nelem();

  out2 << "  Allocates memory for a weighting function matrix having:\n" 
       << "    " << ny << " frequency points\n"
       << "    " << nx << " state variables\n"
       << "    " << ni << " retrieval identities\n";

  kx.resize(         ny, nx );
  kx_names.resize(   ni     );
  kx_lengths.resize( ni     );
  kx_aux.resize(     nx, 2  );

  for ( Index i=0; i<ni; i++ )
    kx_lengths[i] = 0;
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-21
*/
void kbAllocate (
                    Matrix&          kb,
                    ArrayOfString&   kb_names,
                    ArrayOfIndex&    kb_lengths,
                    Matrix&          kb_aux,
              const Vector&          y,
              const String&          y_name,
              const Index&             ni,
              const Index&             nx )
{
  kxAllocate( kb, kb_names, kb_lengths, kb_aux, y, y_name, ni, nx );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-21
*/
void kxPutInK (
		    Matrix&          kx,
		    ArrayOfString&   kx_names,
		    ArrayOfIndex&    kx_lengths,
		    Matrix&          kx_aux,
              const Matrix&          k,
              const ArrayOfString&   k_names,
              const Matrix&          k_aux )
{
  const Index  ny  = kx.nrows();
  const Index  nx  = kx.ncols();
  const Index  ni  = kx_names.nelem();
  const Index  nx2 = k.ncols();
  const Index  ni2 = k_names.nelem();

  if ( ny != k.nrows() )
    throw runtime_error("The two WF matrices have different number of rows.");

  // Determine old length of x and number of identities
  Index  ni1, nx1=0;
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
  Index l = (Index) floor(nx2/ni2);

  // Move K to Kx
  kx(     Range(joker),   Range(nx1,nx2) ) = k;
  kx_aux( Range(nx1,nx2), Range(joker)   ) = k_aux;
  // For the Array types kx_names and kx_lenghts we cannot use Range.  
  for ( Index iri=0; iri<ni2; iri++ )
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
		    Matrix&          kb,
		    ArrayOfString&   kb_names,
		    ArrayOfIndex&    kb_lengths,
		    Matrix&          kb_aux,
              const Matrix&          k,
              const ArrayOfString&   k_names,
              const Matrix&          k_aux )
{
  kxPutInK( kb, kb_names, kb_lengths, kb_aux, k, k_names, k_aux );
}
