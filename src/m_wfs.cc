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

/*-----------------------------------------------------------------------
FILE:      m_wfs.cc

INCLUDES:  

FUNCTIONS: 

HISTORY:   
-----------------------------------------------------------------------*/

#include "arts.h"
#include "vecmat.h"
#include "messages.h"          
#include "workspace.h"          
#include "math_funcs.h"          
#include "atm_funcs.h"          



//==========================================================================
//=== Help functions for grid conversions
//==========================================================================

// GRID2GRID_INDEX
//
// This function gives the indeces of the grid points affected by a change at 
// a set of retrieval points. That is, a first step to calculate dkappa/dkp
// (see Equation 4 of the AUG section on "Basic atmospheric WFs). The second
// stap is performed by the function GRID2GRID_WEIGHTS found below.
//
// This function can be used to convert LOS WFs to WFs for vertical altitudes.
// 
// The profile/function for the retrieved quantity is treated to be piecewise
// linear. With other words, the base functions are tenth functions. The 
// functions areis treated to be 0 outside the end points. 
// The function returns a 2 column matrix, where column 1 is the lowest
// index affected and column 2 the highest index. If no points are affected
// 0's are returned. 
//
// An example:
// x0 = [1 2 3 4 5 6 7] and xp = [0 1 2 3 4.5 5 5.5 6 8] gives
// Is = [ 0     0
//        1     1
//        2     2
//        3     4
//        4     4
//        5     5
//        0     0
//        6     7
//        7     7 ]
// Note that the XP point 5.5 gives 0 as this point not embraces any point 
// of X0 as there are retrieval points at the neighbouring X0 points.
// 
// X0 is the original grid. For example, the LOS grid. 
// XP is the retrieval grid. For example, vertical altitudes
//
// Both vectors most be sorted with increasing values. 
// The length of XP must be > 1.
//
// The algorithm used here is the most effective. However, more elaborated
// algorithms were tested but they failed for some special case, and 
// instead a simplistic approach was used.
// 
//
// Patrick Eriksson 14.06.00

void grid2grid_index (
                    MATRIX&   Is,
              const VECTOR&   x0,
              const VECTOR&   xp )
{
  const size_t n0 = x0.dim();        // length of original grid
  const size_t np = xp.dim();        // length if retrieval grid
        size_t i0, ip;               // counter for each grid

  // Resize Is and set all values to 0
  Is.newsize(np,2);
  Is = 0.0;
 
  // Some checks
  if ( np < 2 )
    throw logic_error("The retrieval grid must have a length > 1.");
  for ( i0=2; i0<=n0; i0++ )
  {
    if ( x0(i0-1) >= x0(i0) )
      throw logic_error("The grids must be sorted with increasing values.");
  }
  for ( ip=2; ip<=np; ip++ )
  {
    if ( xp(ip-1) >= xp(ip) )
      throw logic_error("The grids must be sorted with increasing values.");
  }

  // Do something only if at least one point of X0 is inside XP
  if ( !( (x0(1)>xp(np)) || (x0(n0)<xp(1)) ) )
  {

    // First point of XP
    i0 = 1;
    for ( ; x0(i0) < xp(1); i0++ ) {}
    if ( x0(i0) < xp(2) )
    {
      Is(1,1) = (Numeric) i0;
      for ( ; (i0<=n0) && (x0(i0)<xp(2)); i0++ ) {}
      Is(1,2) = (Numeric) (i0 - 1);
    }

    // Points inside XP
    for ( ip=2; ip<np; ip++ )
    {
      i0 = 1;
      for ( ; (i0<=n0) && (x0(i0)<=xp(ip-1)); i0++ ) {}
      if ( (i0<=n0) && (x0(i0)<xp(ip+1)) )
      {
        Is(ip,1) = (Numeric) i0;
        for ( ; (i0<=n0) && (x0(i0)<xp(ip+1)); i0++ ) {}
        Is(ip,2) = (Numeric) (i0 - 1);
      }
    }

    // Last point of XP
    i0 = 1;
    for ( ; (i0<=n0) && (x0(i0)<=xp(np-1)); i0++ ) {}
    if ( (i0<=n0) && (x0(i0)<=xp(np)) )
    {
      Is(np,1) = (Numeric) i0;
      for ( ; (i0<=n0) && (x0(i0)<=xp(np)); i0++ ) {}
      Is(np,2) = (Numeric) (i0 - 1);
    }
  }
}



// GRID2GRID_INDEX
//
// This function returns the derivative dkappa/dkp (see Equation 4 of 
// the AUG section on "Basic atmospheric WFs"), here dennoted as the weight,
// for the indeces selected.
// This is an accompanying function to GRID2GRID_INDEX. This latter function
// gives the indeces for which the weight > 0.
//
// The weight is
// w = abs( (x0-xp(ip+-1)) / (xp(ip)-xp(ip+-1)));
// where x0 is the point of the original grid (kappa), xp(ip) is the 
// considered retrieval point and xp(ip+-1) is the retrieval point on the
// other side of x0.
//
// X0 is the original grid. For example, the LOS grid. 
// XP is the retrieval grid. For example, vertical altitudes
// I1 is the first index of X0 to consider
// I2 is the last index of X0 to consider
// IP is the index of the retrieval point
//
// Patrick Eriksson 14.06.00

void grid2grid_weights (
                    VECTOR&   w,
              const VECTOR&   x0,
              const size_t&   i1,
              const size_t&   i2,
              const VECTOR&   xp,
              const size_t&   ip )
{
  const size_t   np =xp.dim();        // number of retrieval points
  const size_t   nw = i2-i1+1;        // number of LOS points affected
        size_t   i;

  // Reallocate w
  w.newsize(nw);

  // First point of the retrieval grid
  if ( ip == 1 )
  {
    for ( i=0; i<nw; i++ )
      w(i+1) = ( xp(2)- x0(i1+i) ) / ( xp(2) - xp(1) );
  }
  // Points inside the retrieval grid
  else if ( ip < np )  
  {
    for ( i=0; i<nw; i++ )
    {
      if ( x0(i1+i) <= xp(ip) )
        w(i+1) = ( x0(i1+i) - xp(ip-1) ) / ( xp(ip) - xp(ip-1) );
      else
        w(i+1) = ( xp(ip+1) - x0(i1+i) ) / ( xp(ip+1) - xp(ip) );
    }
  }
  // Last point of the retrieval grid
  else
  {
    for ( i=0; i<nw; i++ )
      w(i+1) = ( x0(i1+i) - xp(np-1) ) / ( xp(np) - xp(np-1) );
  }
}



//==========================================================================
//=== Help functions for KLOS1D
//==========================================================================

// KLOS_1PASS
//
// This function covers cases where each point is passed once.
// That is 1D limb sounding and downward observatiuons are not covered.
//
// The expression used are described in sub-section 2.1 of the AUG section
// "Basic atmospheric WFs".
//
// Patrick Eriksson 09.06.00

void klos_1pass (
                    MATRIX&   K,
                    VECTOR    y,
              const int&      start_index,
              const int&      stop_index,
              const Numeric&  lstep,
              const MATRIX&   Tr,
              const MATRIX&   S,
              const int&      ground,
              const VECTOR&   e_ground,
              const VECTOR&   y_ground )
{
  const size_t   nf = Tr.dim(1);      // number of frequencies
        size_t   iv;                  // frequency index
        VECTOR   t(nf,1.0);           // transmission to the sensor
           int   q;                   // index corresponding to point q in AUG

  if ( (ground==1) || (ground==start_index) )
    throw logic_error("The ground cannot be at one of the end points."); 

  // Resize K
  K.newsize( nf, start_index );

  // We start here at the LOS point closest to the sensor, that is,
  // reversed order compared to RTE_ITERATE  

  // The LOS point closest to the sensor 
  q  = stop_index;
  for ( iv=1; iv<=nf; iv++ )    
  {
    t(iv)   *= Tr(iv,q);
    y(iv)    = (y(iv)-S(iv,q)*(1.0-Tr(iv,q)))/Tr(iv,q);
    K(iv,q)  = -lstep*(y(iv)-S(iv,q))*t(iv)/2;
  }

  // Points inside the LOS
  for ( q=stop_index+1; q<start_index; q++ )
  {
    for ( iv=1; iv<=nf; iv++ )    
    {
      y(iv)   = (y(iv)-S(iv,q)*(1.0-Tr(iv,q)))/Tr(iv,q);
      K(iv,q) = -lstep*(2*(y(iv)-S(iv,q))*Tr(iv,q)+S(iv,q)-S(iv,q-1))*t(iv)/2;
      t(iv)  *= Tr(iv,q); 
    }
  }

  // The LOS point furthest away from the sensor
  q = start_index;
  for ( iv=1; iv<=nf; iv++ )    
    K(iv,q)  = -lstep*(y(iv)-S(iv,q-1))*t(iv)/2;

  // Y shall now be equal to the radiation entering the atmosphere and T
  // the total transmission
}



// KLOS_1DLIMB
//
// This function covers 1D limb sounding
//
// The expression used are described in sub-section 2.2 of the AUG section
// "Basic atmospheric WFs".
//
// Patrick Eriksson 14.06.00

void klos_1dlimb (
                    MATRIX&   K,
                    VECTOR    y,
                    VECTOR    yn,             // = y_space
              const int&      start_index,
              const Numeric&  lstep,
              const MATRIX&   Tr,
              const MATRIX&   S,
              const int&      ground,
              const VECTOR&   e_ground )
{
  const size_t   nf = Tr.dim(1);      // number of frequencies
        size_t   iv;                  // frequency index
        VECTOR   t1q;                 // transmission tangent point - q
        VECTOR   tqn(nf,1);           // transmission q - sensor
        int      q;                   // index corresponding to point q in AUG
        Numeric  tv, tv1;             // transmission value for q and q-1

  // Calculate the square root of the total transmission
  bl( t1q, start_index, start_index, Tr, ground, e_ground );
  t1q = sqrt(t1q);

  // Resize K
  K.newsize( nf, start_index );

  // We start at the outermost point
  q  = start_index;       
  for ( iv=1; iv<=nf; iv++ )    
  {
    tv1      = Tr(iv,q-1);
    t1q(iv) /= tv1;
    tqn(iv) *= tv1;
    y(iv)    = ( y(iv) - S(iv,q-1)*(1-tv1)*(1+t1q(iv)*t1q(iv)*tv1) ) / tv1;
    K(iv,q)  = -lstep*( ( 2*yn(iv) + S(iv,q-1)*(1-2*tv1) )*t1q(iv)*t1q(iv)*tv1 
                        + y(iv) - S(iv,q-1) )*tv1/2;
  }

  // Points inside the LOS
  for ( q=start_index-1; q>1; q-- )
  {
    for ( iv=1; iv<=nf; iv++ )    
    {
      tv1      = Tr(iv,q-1);    
      tv       = Tr(iv,q);
      t1q(iv) /= tv1;
      y(iv)    = ( y(iv) - S(iv,q-1)*(1-tv1)*(1+t1q(iv)*t1q(iv)*tv1) ) / tv1;
      K(iv,q)  = -lstep*( ( 4*(yn(iv)-S(iv,q)) + 3*(S(iv,q)-S(iv,q-1)) + 
                            2*S(iv,q-1) )*t1q(iv)*t1q(iv)*tv1 + 
                     2*(y(iv)-S(iv,q-1))*tv1 + S(iv,q-1) - S(iv,q) )*tqn(iv)/2;
      tqn(iv) *= tv1;
      yn(iv)   = yn(iv)*tv + S(iv,q)*(1-tv);
    } 
  }

  // The tangent or ground point
  for ( iv=1; iv<=nf; iv++ )    
    K(iv,q)  = -lstep*( y(iv) - S(iv,1) )*tqn(iv)*Tr(iv,1);

}



// KLOS_1DDOWN
//
// This function covers 1D downward looking observations
//
// The expression used are described in sub-section 2.3 of the AUG section
// "Basic atmospheric WFs".
//
// Patrick Eriksson 15.06.00

void klos_1ddown (
                    MATRIX&   K,
              const VECTOR    y,
              const VECTOR    y_space,
              const int&      start_index,
              const int&      stop_index,
              const Numeric&  lstep,
              const MATRIX&   Tr,
              const MATRIX&   S,
              const int&      ground,
              const VECTOR&   e_ground,
              const VECTOR&   y_ground )
{
  const size_t   nf = Tr.dim(1); // number of frequencies
        size_t   iv;             // frequency index
           int   q;              // LOS index
        VECTOR   y0;             // see below
        MATRIX   K2;             // matrix for calling other LOS WFs functions
        VECTOR   tr0;            // see below

  // Resize K
  K.newsize( nf, start_index );

  // Calculate Y0, the intensity reaching the platform altitude at the far
  // end of LOS, that is, from above.
  y0 = y_space;
  rte_iterate( y0, start_index-1, stop_index, Tr, S, nf );

  // Calculate TR0, the transmission from the platform altitude down to the
  // tangent point or the ground, and up to the platform again.
  bl( tr0, stop_index, stop_index, Tr, ground, e_ground );


  // The indeces below STOP_INDEX are handled by the limb sounding function.
  // The limb function is given Y0 instead of cosmic radiation 
  klos_1dlimb( K2, y, y0, stop_index, lstep, Tr, S, ground, e_ground );  
  for ( iv=1; iv<=nf; iv++ )
  {
    for ( q=1; q<stop_index; q++ )
      K(iv,q) = K2(iv,q);
  }

  // The indeces above STOP_INDEX are handled by the 1pass function.
  // The transmission below STOP_INDEX must here be considered.
  klos_1pass( K2, y0, start_index, stop_index, lstep, Tr, S, 
                     ground, e_ground, y_ground );
  for ( iv=1; iv<=nf; iv++ )
  {
    for ( q=stop_index+1; q<=start_index; q++ )
      K(iv,q) = K2(iv,q)*tr0(iv);
  }

  // The platform altitude must be treated seperately
  //
  // Calculate the intensity generated below point q-1, YQQ
  VECTOR yqq;
  rte( yqq, stop_index-1, stop_index-1, Tr, S, VECTOR(nf,0.0), 
                                            ground, e_ground, y_ground );
  //
  // Y0 is moved one step upwards and TR0 one step downwards
  q = stop_index; 
  for ( iv=1; iv<=nf; iv++ )
  {
    tr0(iv) /= Tr(iv,q-1)*Tr(iv,q-1);    
    y0(iv)   = ( y0(iv) - S(iv,q)*(1-Tr(iv,q)) ) / Tr(iv,q);
    K(iv,q)  = -lstep*( (3*(y0(iv)-S(iv,q))*Tr(iv,q-1)*Tr(iv,q) + 
                      2*(S(iv,q)-S(iv,q-1))*Tr(iv,q-1) + S(iv,q-1) )*tr0(iv) + 
                                           yqq(iv) - S(iv,q-1) )*Tr(iv,q-1)/2;
  }

}

//==========================================================================
//=== Workspace methods
//==========================================================================

// KLOS1D
//
// Patrick Eriksson 14.06.00

void klos1d (
                    ARRAYofMATRIX&   k,
              const Los&             los,   
              const ARRAYofMATRIX&   source,
              const ARRAYofMATRIX&   trans,
              const VECTOR&          y,
              const VECTOR&          y_space,
              const VECTOR&          f_abs,
              const VECTOR&          e_ground,
              const Numeric&         t_ground )
{
  const size_t  nfi = los.start.dim();   // number of viewing angles  
  const size_t  nf  = f_abs.dim();       // number of frequencies  
        VECTOR  yp(nf);                  // part of Y
        size_t  iy, iy0=0;               // y indices

  // Set up vector for ground blackbody radiation
  VECTOR   y_ground; 
  if ( any(los.ground) )
    planck( y_ground, f_abs, t_ground );

  // Resize the LOS WFs array
  k.newsize(nfi);

  // Loop viewing angles
  for (size_t i=1; i<=nfi; i++ ) 
  {
    // Do something only if LOS has any points
    if ( los.p(i).dim() > 0 )
    {
      // Pick out the part of Y corresponding to the present viewing angle
      for ( iy=1; iy<=nf; iy++ )    
        yp(iy) = y(iy0+iy);
      iy0 += nf;        

      // The calculations are performed in 3 sub-functions
      //
      // Each point only passed once
      if ( los.stop(i)==1 )
        klos_1pass( k(i), yp, los.start(i), 1, los.l_step(i), 
                     trans(i), source(i), los.ground(i), e_ground, y_ground );
      //
      // 1D limb sounding
      else if ( los.start(i) == los.stop(i) )
        klos_1dlimb( k(i), yp, y_space, los.start(i), los.l_step(i), 
                     trans(i), source(i), los.ground(i), e_ground );
      //
      // 1D downward looking
      else 
        klos_1ddown( k(i), yp, y_space, los.start(i), los.stop(i), 
                        los.l_step(i), trans(i), source(i), los.ground(i), 
                        e_ground, y_ground );
    }
  }
}



// KSPECIES1D
//
// Patrick Eriksson 14.06.00

void kSpecies1d (
                    MATRIX&          K,
              const Los&             los,           
              const ARRAYofMATRIX&   klos,
              const VECTOR&          p_abs,
              const MATRIX&          Abs,
              const VECTOR&          p_grid )
{
  const size_t  nfi = los.start.dim();   // number of viewing angles  
  const size_t  nf  = Abs.dim(1);        // number of frequencies
  const size_t  np  = p_grid.dim();      // number of retrieval points  
        MATRIX  Abs2;                    // absorption at p_grid
        VECTOR  lplos;                   // -log of the LOS pressures  
  const VECTOR  lgrid = -1.0*log(p_grid);// -log of the grid pressures
        MATRIX  Is;                      // matrix for storing LOS index
        VECTOR  w;                       // weights for LOS WFs
        VECTOR  a(nf);                   // temporary vector
        size_t  ip;                      // Retrieval point index
        size_t  iv, iv0=0;               // Frequency index
        size_t  i1, iw;                  // weight index

  // Resize K and set all values to 0
  K.newsize(nfi*nf,np);
  K = 0.0;

  // Determine the absorption at the retrieval points
  interpp( Abs2, p_abs, Abs, p_grid );

  // Loop viewing angles
  for (size_t i=1; i<=nfi; i++ ) 
  {
    // Do something only if LOS has any points
    if ( los.p(i).dim() > 0 )
    {
      // Get the LOS points affected by each retrieval point
      lplos = -1.0 * log(los.p(i));
      grid2grid_index( Is, lplos, lgrid );

      // Loop retrieval points
      for ( ip=1; ip<=np; ip++ ) 
      {
        // Check if there is anything to do
        if ( Is(ip,1) > 0 )
        {
          // Get the weights for the LOS points
	  grid2grid_weights( w, lplos, (size_t)Is(ip,1), (size_t) Is(ip,2), 
                                                                 lgrid, ip );
          // Calculate the WFs.
          // A is di/dkappa*dkappa/dkp in a compact form.
          // This is possible as the columns of dkappa/dkp are identical.  
          a = 0.0;                     
          i1 = (size_t)Is(ip,1);       // first LOS point to consider
          for ( iv=1; iv<=nf; iv++ )
	  {
            for ( iw=i1; iw<=(size_t)Is(ip,2); iw++ )
              a(iv) += klos(i)(iv,iw) * w(iw-i1+1);
            K(iv0+iv,ip) = a(iv) * Abs2(iv,ip);                    
	  }
        }            
      }
    }
     
    // Update the frequency index offset
    iv0 += nf;
  }  
}
