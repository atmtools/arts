/* Copyright (C) 2002 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                 2002 Oliver Lemke <olemke@sat.physik.uni-bremen.de>
     
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
   USA. 
*/

/*!
  \file   gridded_fields.cc
  \author Oliver Lemke  <olemke@sat.physik.uni-bremen.de>
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
  \date   Wed Jun 26 17:48:29 2002
  
  \brief  Reading routines for gridded fields.
  
  This file contains reading routines for gridded fields. Gridded fields are
  needed to store moredimesional data together with the corresponding grids 
  in the same variable. The datatype og a gridded field is always an array
  of tensors.

  For further description see AUG.

*/

/*===========================================================================
  === External declarations
  ===========================================================================*/


#include <stdexcept>
#include <iostream>
#include "xml_io.h"
#include "array.h"
#include "matpackIII.h"
#include "matpackVI.h"
#include "gridded_fields.h"

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

ostream& operator<< (ostream &os, const GriddedField3 &/*gfield3*/)
{
  os << "GriddedField3: Output operator not implemented";
  return os;
}


ostream& operator<< (ostream &os, const ArrayOfGriddedField3 &/*agfield3*/)
{
  os << "ArrayOfGriddedField3: Output operator not implemented";
  return os;
}

//! Check the dimensions of a gridded tensor3
/*! 
  
  \param gridded_tensor Gridded tensor3.
*/
void
check_gridded_tensor3 (ArrayOfTensor3 &gridded_tensor)
{
  // If the field has 3 dimensions the gridded field tensor must have 4
  // elements. For each dimension there is 1 grid. All grids are stored in
  // a Tensor3. So there are 3 Tensor3 for the grids and additionally the
  // Tensor3 containing the data.
  if(gridded_tensor.nelem() != 4)
    throw runtime_error("Error while reading data file (gridded field format):"
                        "A gridded Tensor3 must have 4 elements.");

  
  // Create array of functions to be able to retrieve the tensor dimensions
  // index instead of name. 0 = nvitrines(), 1 = nshelves(), 2 = nbooks() ...
  Index (Tensor3::*f[3])() const =
    {
      &Tensor3::npages,
      &Tensor3::nrows,
      &Tensor3::ncols
    };

  // Create array with the dimension names, so we can access them
  // easily by index similar to the tensor dimensions before. Used for
  // error output.
  static const char elem_names[][6] =
    { 
      "pages",
      "rows",
      "cols" };
  
  // Loop over all grid tensors
  for (Index g=0; g < 3; g++)
    {
      bool grid_empty = true;

      // Check whether the grid tensor is empty (all dimensions are 0)
      for (Index d=0; grid_empty && d < 3; d++)
        {
          if ((gridded_tensor[g].*f[d])() != 0)
            grid_empty = false;
        }

      if (!grid_empty)
        {
          // Loop over all grid tensor dimensions
          for (Index d=0; d < 3; d++)
            {
              if (g == d)
                {
                  // If the grid tensor is the same as the grid tensor
                  // dimension, the corresponding dimension of the
                  // data tensor must be the same
                  if ((gridded_tensor[g].*f[d])() !=
                      (gridded_tensor[3].*f[d])())
                    {
                      ostringstream os;
                      os <<"Error while reading data file "
                         << "(gridded data field format): \n"
                         <<elem_names[d] << " dimension of "
                         << elem_names[d] << "_grid "
                         << "(" << (gridded_tensor[g].*f[d])() << ") "
                         << "must be the same as in data tensor "
                         << "(" << (gridded_tensor[3].*f[d])() << ") "
                         << "\n";
                      throw runtime_error (os.str());
                    }
                }
              else
                {
                  // All other dimensions of the grid tensor must be 1
                  if ((gridded_tensor[g].*f[d])() != 1)
                    {
                      ostringstream os;
                      os << "Error while reading data file "
                         << "(gridded data field format): \n"
                         << elem_names[d] << " dimension of "
                         << elem_names[g] << "_grid "
                         << "(" << (gridded_tensor[g].*f[d])() << ") "
                         << "must be 1"
                         << "\n";
                      throw runtime_error (os.str());
                    }
                }
            }
        }
      else
        {
          // If the grid tensor is empty the corresponding dimension of
          // the data tensor must be 1
          if ((gridded_tensor[3].*f[g])() != 1)
            {
              ostringstream os;
              os << "Error while reading data file "
                 << "(gridded data field format): \n"
                 << "If the " << elem_names[3] << "_grid is empty, the "
                 << elem_names[3] << " dimension "
                 << "(" << (gridded_tensor[3].*f[g])() << ") "
                 << "of the data tensor must be 1"
                 << "\n";
              throw runtime_error (os.str());
            }
        }
    }

}

//! Check the dimensions of a gridded tensor6
/*! 
  
  \param gridded_tensor Gridded tensor6.
*/
void
check_gridded_tensor6 (ArrayOfTensor6 &gridded_tensor)
{

  // If the field has 3 dimensions the gridded field tensor must have 4
  // elements. For each dimension there is 1 grid. All grids are stored in
  // a Tensor3. So there are 3 Tensor3 for the grids and additionally the
  // Tensor3 containing the data.
  if(gridded_tensor.nelem() != 7)
    throw runtime_error("Error while reading data file (gridded field format):"
                        "A gridded Tensor3 must have 4 elements."); 

  // Create array of functions to be able to retrieve the tensor dimensions
  // index instead of name. 0 = nvitrines(), 1 = nshelves(), 2 = nbooks() ...
  Index (Tensor6::*f[6])() const =
    {
      &Tensor6::nvitrines,
      &Tensor6::nshelves,
      &Tensor6::nbooks,
      &Tensor6::npages,
      &Tensor6::nrows,
      &Tensor6::ncols
    };

  // Create array with the dimension names, so we can access them
  // easily by index similar to the tensor dimensions before. Used for
  // error output.
  static const char elem_names[][9] =
    { "vitrines",
      "shelves",
      "books",
      "pages",
      "rows",
      "cols" };
  
  // Loop over all grid tensors
  for (Index g=0; g < 6; g++)
    {
      bool grid_empty = true;

      // Check whether the grid tensor is empty (all dimensions are 0)
      for (Index d=0; grid_empty && d < 6; d++)
        {
          if ((gridded_tensor[g].*f[d])() != 0)
            grid_empty = false;
        }

      if (!grid_empty)
        {
          // Loop over all grid tensor dimensions
          for (Index d=0; d < 6; d++)
            {
              if (g == d)
                {
                  // If the grid tensor is the same as the grid tensor
                  // dimension, the corresponding dimension of the
                  // data tensor must be the same
                  if ((gridded_tensor[g].*f[d])() !=
                      (gridded_tensor[6].*f[d])())
                    {
                      ostringstream os;
                      os << "Error while reading data file "
                         << "(gridded data field format): \n"
                         << elem_names[d] << " dimension of "
                         << elem_names[d] << "_grid "
                         << "(" << (gridded_tensor[g].*f[d])() << ") "
                         << "must be the same as in data tensor "
                         << "(" << (gridded_tensor[6].*f[d])() << ") "
                         << "\n";
                      throw runtime_error (os.str());
                    }
                }
              else
                {
                  // All other dimensions of the grid tensor must be 1
                  if ((gridded_tensor[g].*f[d])() != 1)
                    {
                      ostringstream os;
                      os << "Error while reading data file "
                         << "(gridded data field format): \n"
                         << elem_names[d] << " dimension of "
                         << elem_names[g] << "_grid "
                         << "(" << (gridded_tensor[g].*f[d])() << ") "
                         << "must be 1"
                         << "\n";
                      throw runtime_error (os.str());
                    }
                }
            }
        }
      else
        {
          // If the grid tensor is empty the corresponding dimension of
          // the data tensor must be 1
          if ((gridded_tensor[6].*f[g])() != 1)
            {
              ostringstream os;
              os << "Error while reading data file "
                 << "(gridded data field format): \n"
                 << "If the " << elem_names[g] << "_grid is empty, the "
                 << elem_names[g] << " dimension "
                 << "(" << (gridded_tensor[6].*f[g])() << ") "
                 << "of the data tensor must be 1"
                 << "\n";
              throw runtime_error (os.str());
            }
        }
    }

}



//! This function reads a 3D gridded field. 
/*! 
  A 3D gridded field is of datatype ArrayOfTensor3. The array consists of 4 
  elements. The first three elements contain the grids and the fourth the data.
  The dimensions of the four tensors are:

  'pages'-grid:  [N_pages_grid, 1, 1]
  'rows'-grid:   [1, N_rows_grid, 1]
  'cols'-grid:   [1, 1, N_cols_grid]
  data:          [N_pages_grids, N_rows_grid, N_cols_grid]
  
  where N_pages_grid, N_row_grid and N_cols_grid are the sizes of the grids 
  on which the data is stored on.

  The data files must be in xml-format. This function reads a gridded field
  of dimension 3. It uses the xml reading routine for an ArrayOfTensor3 and 
  checks if the data is in a 'gridded field' format as described above.

  \param filename         Name of the file to be read.
  \param gridded_tensor3  Name of the variable used for storing the data.
  
 */
void read_gridded_tensor3(
                          const String& filename,
                          ArrayOfTensor3& gridded_tensor3
                          )
{ 
  // Read the griddes field from a file.
  xml_read_from_file( filename, gridded_tensor3 );

  // Check if the tensor is a "gridded" tensor.
  check_gridded_tensor3(gridded_tensor3);
}


//! This function reads a 6D gridded field. 
/*! 
  A 6D gridded field is of datatype ArrayOfTensor6. The array consists of 7 
  elements. The first six elements contain the grids and the seventhth the
  data.
  The dimensions of the seven tensors are:
  'vitrines'-grid:  [N_vitrines_grid, 1, 1, 1, 1, 1]
  'shelves'-grid:   [1, N_shelves_grid, 1, 1, 1, 1]
  'books'-grid:     [1, 1, N_books_grid, 1, 1, 1]
  'pages'-grid:     [1, 1, 1, N_pages_grid, 1, 1]
  'rows'-grid:      [1, 1, 1, 1, N_rows_grid, 1]
  'cols'-grid:      [1, 1, 1, 1, 1, N_cols_grid]
  data:             [N_pages_grids, N_rows_grid, N_cols_grid]
  
  where N_pages_grid, N_row_grid and N_cols_grid etc. are the sizes of the 
  grids on which the data is stored on.

  The data files must be in xml-format. This function reads a gridded field
  of dimension 3. It uses the xml reading routine for an ArrayOfTensor6 and 
  checks if the data is in a 'gridded field' format as described above using 
  the *check_gridded_tensor6* function.

  \param filename         Name of the file to be read.
  \param gridded_tensor6  Name of the variable used for storing the data.
  
 */ 
void read_gridded_tensor6(
                          const String& filename,
                          ArrayOfTensor6& gridded_tensor6
                          )
{
  // Read the griddes field from a file.
  xml_read_from_file( filename, gridded_tensor6 );
  
  // Check if the tensor is a "gridded" tensor.
  check_gridded_tensor6(gridded_tensor6);
}
 
   
