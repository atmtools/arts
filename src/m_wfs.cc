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
  const int      nf = Tr.dim(1);      // number of frequencies
  const int      np = Tr.dim(2)+1;    // number of LOS points
        int      iv;                  // frequency index
        int      q = stop_index;      // LOS point index
        int      q1;                  // one step offset from q
        int      io;                  // index off-set, -1 or 0
        int      is;                  // step order in loop , -1 or 1
        VECTOR   t(nf,1.0);           // transmission to the sensor

  // Resize K
  K.newsize( nf, np );

  // Check loop order
  // We start here at the LOS point closest to the sensor, that is,
  // reversed order compared to RTE_ITERATE  
  if ( start_index == 1 )
  {
    io  = -1;     // We start here at the highest index
    is  = -1;
  }
  else
  {
    io  = 0;      // We start here at index 1
    is  = 1;
  } 

  // Do the calculations
  // The LOS point closest to the sensor (q already set)
  q1 = q+io;
  for ( iv=1; iv<=nf; iv++ )    
  {
    t(iv)   *= Tr(iv,q1);
    y(iv)    = (y(iv)-S(iv,q1)*(1.0-Tr(iv,q1)))/Tr(iv,q1);
    K(iv,q)  = -lstep*(y(iv)-S(iv,q1))*t(iv)/2;
  }
  //
  // Points inside the LOS
  for ( q=stop_index+is; q!=start_index; q+=is )
  {
    q1 = q+io;
    for ( iv=1; iv<=nf; iv++ )    
    {
      K(iv,q) = -lstep*(2*(y(iv)-S(iv,q1))*Tr(iv,q1)+S(iv,q1)-S(iv,q))*t(iv)/2;
      t(iv)  *= Tr(iv,q1); 
      y(iv)   = (y(iv)-S(iv,q1)*(1.0-Tr(iv,q1)))/Tr(iv,q1);
      // Check for ground
    }
  }
  //
  // The LOS point furthest away from the sensor (q1 has correct value)
  for ( iv=1; iv<=nf; iv++ )    
    K(iv,q)  = -lstep*(y(iv)-S(iv,q1))*t(iv)/2;
}



void klosBasic (
                    ARRAYofMATRIX&   k,
              const Los&             los,   
              const ARRAYofMATRIX&   source,
              const ARRAYofMATRIX&   trans,
              const VECTOR&          y,
              const VECTOR&          f_abs,
              const VECTOR&          e_ground,
              const Numeric&         t_ground )
{
  const size_t  nfi = los.start.dim();   // number of viewing angles  
  const size_t  nf  = los.start.dim();   // number of frequencies  
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
    if ( los.p(i).dim() > 0 )
    {
      // Pick out the part of Y corresponding to the present viewing angle
      for ( iy=1; iy<=nf; iy++ )    
        yp(iy) = y(iy0+iy);
      iy0 += nf;        

      // The calculations are performed in 3 sub-functions
      //
      // Each point only passed once
      if (  (los.start(i)==1) || (los.stop(i)==1) )
        klos_1pass( k(i), yp, los.start(i), los.stop(i), los.l_step(i), 
                     trans(i), source(i), los.ground(i), e_ground, y_ground );
      //
      /*      // 1D limb sounding
      else if ( los.start(i) == los.stop(i) )
        klos_1Dlimb(...);
      //
      // 1D downward looking
      else 
      klos_1Dup(...); */
    }
  }
}


/*
void kSpecies (
                    MATRIX&          k,
              const Los&             los,           
              const ARRAYofMATRIX&   klos,
              const VECTOR&          p_abs,
              const MATRIX&          abs,
              const VECTOR&          p_grid )
{
  const size_t  nfi = los.start.dim();   // number of viewing angles  
  const size_t  nf  = abs.dim(1);        // number of frequencies
  const size_t  np  = p_grid.dim();      // number of retrieval points  
        VECTOR  abs2;                    // absorption at p_grid

  // Resize K and set all values to 0
  k.newsize(nfi*nf,np);
  k = 0;

  // Determine the absorption at the retrieval points
  interp_lin_row( abs2, p_abs, abs, p_grid );

  // Loop viewing angles
  for (size_t i=1; i<=nfi; i++ ) 
  {

  }  
}
*/
