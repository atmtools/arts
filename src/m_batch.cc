/* Copyright (C) 2004 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>
                            
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



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_batch.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2004-09-15 

  \brief  Workspace functions for doing batch calculations.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts.h"
#include "auto_md.h"
#include "math_funcs.h"



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! BatchUpdateMatrix
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2004-09-15
*/
void BatchUpdateMatrix(
      // WS Generic Output:
      Matrix&          m,
      // WS Generic Output Names:
      const String&    m_name,
      // WS Input:
      const Index&     ybatch_index,
      // WS Generic Input:
      const Tensor3&   t3,
      // WS Generic Input Names:
      const String&    t3_name )
{
  if( ybatch_index >= t3.npages() )
    {
      ostringstream os;
      os << "The value of *ybatch_index* (" << ybatch_index 
         << "is outside the page range of *" << t3_name << "*.";
      throw runtime_error( os.str() );

    }
  out3 << "Copies page " << ybatch_index << " of *" << t3_name
       << "* to create *" << m_name << "*.\n";

  m.resize( t3.nrows(), t3.ncols() );
  m = t3( ybatch_index, joker, joker );
}



//! BatchUpdateNumeric
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2004-09-15
*/
void BatchUpdateNumeric(
      // WS Generic Output:
      Numeric&         n,
      // WS Generic Output Names:
      const String&    n_name,
      // WS Input:
      const Index&     ybatch_index,
      // WS Generic Input:
      const Vector&    v,
      // WS Generic Input Names:
      const String&    v_name )
{
  if( ybatch_index >= v.nelem() )
    {
      ostringstream os;
      os << "The value of *ybatch_index* (" << ybatch_index 
         << "is outside the range of the vector *" << v_name << "*.";
      throw runtime_error( os.str() );

    }
  out3 << "Copies column " << ybatch_index << " of *" << v_name
       << "* to create *" << n_name << "*.\n";

  n = v[ ybatch_index ];
}



//! BatchUpdateTensor3
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2004-09-15
*/
void BatchUpdateTensor3(
      // WS Generic Output:
      Tensor3&         t3,
      // WS Generic Output Names:
      const String&    t3_name,
      // WS Input:
      const Index&     ybatch_index,
      // WS Generic Input:
      const Tensor4&   t4,
      // WS Generic Input Names:
      const String&    t4_name )
{
  if( ybatch_index >= t4.nbooks() )
    {
      ostringstream os;
      os << "The value of *ybatch_index* (" << ybatch_index 
         << "is outside the book range of *" << t4_name << "*.";
      throw runtime_error( os.str() );

    }
  out3 << "Copies book " << ybatch_index << " of *" << t4_name
       << "* to create *" << t3_name << "*.\n";

  t3.resize( t4.npages(), t4.nrows(), t4.ncols() );
  t3 = t4( ybatch_index, joker, joker, joker );
}



//! BatchUpdateTensor4
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2004-09-15
*/
void BatchUpdateTensor4(
      // WS Generic Output:
      Tensor4&         t4,
      // WS Generic Output Names:
      const String&    t4_name,
      // WS Input:
      const Index&     ybatch_index,
      // WS Generic Input:
      const Tensor5&   t5,
      // WS Generic Input Names:
      const String&    t5_name )
{
  if( ybatch_index >= t5.nshelves() )
    {
      ostringstream os;
      os << "The value of *ybatch_index* (" << ybatch_index 
         << "is outside the shelf range of *" << t5_name << "*.";
      throw runtime_error( os.str() );

    }
  out3 << "Copies shelf " << ybatch_index << " of *" << t5_name
       << "* to create *" << t4_name << "*.\n";

  t4.resize( t5.nbooks(), t5.npages(), t5.nrows(), t5.ncols() );
  t4 = t5( ybatch_index, joker, joker, joker, joker );
}



//! BatchUpdateVector
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2004-09-15
*/
void BatchUpdateVector(
      // WS Generic Output:
      Vector&          v,
      // WS Generic Output Names:
      const String&    v_name,
      // WS Input:
      const Index&     ybatch_index,
      // WS Generic Input:
      const Matrix&    m,
      // WS Generic Input Names:
      const String&    m_name )
{
  if( ybatch_index >= m.nrows() )
    {
      ostringstream os;
      os << "The value of *ybatch_index* (" << ybatch_index 
         << "is outside the row range of *" << m_name << "*.";
      throw runtime_error( os.str() );

    }
  out3 << "Copies row " << ybatch_index << " of *" << m_name
       << "* to create *" << v_name << "*.\n";

  v.resize( m.ncols() );
  v = m( ybatch_index, joker );
}



//! ybatchCalc
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2004-09-15
*/
void ybatchCalc(
        // WS Output:
              Matrix&         ybatch,
              Index&          ybatch_index,
              Index&          ybatch_n,
              Vector&         y,
        // WS Input:
        const Agenda&         batch_update_agenda,
        const Agenda&         batch_calc_agenda,
        const Agenda&         batch_post_agenda )
{
  for( ybatch_index=0; ybatch_index<ybatch_n; ybatch_index++ )
    {
      batch_update_agenda.execute( 0 );
      batch_calc_agenda.execute( 0 );
      
      if( ybatch_index == 0 )
        { ybatch.resize( y.nelem(), ybatch_n); }

      ybatch( joker, ybatch_index ) = y;
    }

  batch_post_agenda.execute( 0 );
}
