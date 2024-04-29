/*--------------------------------------------------------------------
 * $Id: locate.c 2623 2011-12-23 10:52:38Z robert.buras $
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

int locate (double *xx, int n, double x)
/* Given an array xx[1..n], and given a value x, returns a value j such that x is between xx[j] */
/* and xx[j+1]. xx must be monotonic, either increasing or decreasing. j=0 or j=n is returned   */
/* to indicate that x is out of range. A unit-offset array xx is assumed. To use locate with    */
/* a zero-offset array, remember to subtract 1 from the address of xx, and also                 */
/* from the returned value j.                                                                   */
/*                                                                                              */
/* From "Numerical recipes in C"                                                                */
{
  int ju=0, jm=0, jl=0, j=0;
  int monotonicity=0;
  
  jl=0;   /* Initialize lower  */
  ju=n+1; /* and upper limits. */
  
  if (n==1) /* n==1 case is special case */
    monotonicity=1;
  else if (xx[1]>xx[0])
    monotonicity=1;
  
  xx-=1;

  while (ju-jl > 1) { /* If we are not yet done,             */
    jm=(ju+jl) >> 1;  /* compute a midpoint,                 */
    if ((x >= xx[jm] && monotonicity==1) || (x <= xx[jm] && monotonicity==0))
      jl=jm;          /* and replace either the lower limit  */
    else
      ju=jm;          /* or the upper limit, as appropriate. */
  }         /* Repeat until the test condition is satisfied. */
  
  if (x == xx[1]) { 
    j=1; /* Then set the output */
  }
  else {
    if (x == xx[n])
      j=n-1;
    else
      j=jl;
  }

  j-=1;

  return j;
} 


int flocate (float *xx, int n, float x)
/* Given an array xx[1..n], and given a value x, returns a value j such that x is between xx[j] */
/* and xx[j+1]. xx must be monotonic, either increasing or decreasing. j=0 or j=n is returned   */
/* to indicate that x is out of range. A unit-offset array xx is assumed. To use locate with    */
/* a zero-offset array, remember to subtract 1 from the address of xx, and also                 */
/* from the returned value j.                                                                   */
/*                                                                                              */
/* From "Numerical recipes in C"                                                                */
{
  int ju=0, jm=0, jl=0, j=0;
  int monotonicity=0;
  
  jl=0;   /* Initialize lower  */
  ju=n+1; /* and upper limits. */

  if (n==1) /* n==1 case is special case */
    monotonicity=1;
  else if (xx[1]>xx[0])
    monotonicity=1;
  
  xx-=1;

  while (ju-jl > 1) { /* If we are not yet done,             */
    jm=(ju+jl) >> 1;  /* compute a midpoint,                 */
/*    fprintf(stderr,"x %e >= xx %e; ju %d jm %d jl %d\n",x,xx[jm],ju,jm,jl);*/
    if ((x >= xx[jm] && monotonicity==1) || (x <= xx[jm] && monotonicity==0))
      jl=jm;          /* and replace either the lower limit  */
    else
      ju=jm;          /* or the upper limit, as appropriate. */
  }         /* Repeat until the test condition is satisfied. */
  
  if (x == xx[1]) { 
    j=1; /* Then set the output */
  }
  else {
    if (x == xx[n])
      j=n-1;
    else
      j=jl;
  }

  j-=1;

  return j;
} 
