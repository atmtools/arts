#include <iostream>
#include "array.h"
#include "matpackVI.h"

void
check_gridded_tensor6 (ArrayOfTensor6 &gridded_tensor);

int
main (int argc, char *argv[])
{
  ArrayOfTensor6 gridded_tensor;

  gridded_tensor.resize (7);

  gridded_tensor[0].resize (5, 1, 1, 1, 1, 1);
  gridded_tensor[1].resize (1, 4, 1, 1, 1, 1);
  gridded_tensor[2].resize (1, 1, 3, 1, 1, 1);
  gridded_tensor[3].resize (1, 1, 1, 6, 1, 1);
  gridded_tensor[4].resize (1, 3, 1, 1, 6, 1);
  gridded_tensor[5].resize (0, 0, 0, 0, 0, 0);
  gridded_tensor[6].resize (5, 4, 3, 6, 7, 2);

  check_gridded_tensor6 (gridded_tensor);

  return (0);
}

class test
{
public:
  Index nvitrines () const { return 1; }
};

void
check_gridded_tensor6 (ArrayOfTensor6 &gridded_tensor)
{
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
                      cout << "Error while reading data file "
                           << "(gridded data field format): \n"
                           << elem_names[d] << " dimension of "
                           << elem_names[d] << "_grid "
                           << "(" << (gridded_tensor[g].*f[d])() << ") "
                           << "must be the same as in data tensor "
                           << "(" << (gridded_tensor[6].*f[d])() << ") "
                           << endl << endl;
                    }
                }
              else
                {
                  // All other dimensions of the grid tensor must be 1
                  if ((gridded_tensor[g].*f[d])() != 1)
                    {
                      cout << "Error while reading data file "
                           << "(gridded data field format): \n"
                           << elem_names[d] << " dimension of "
                           << elem_names[g] << "_grid "
                           << "(" << (gridded_tensor[g].*f[d])() << ") "
                           << "must be 1"
                           << endl << endl;
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
              cout << "Error while reading data file "
                   << "(gridded data field format): \n"
                   << "If the " << elem_names[g] << "_grid is empty, the "
                   << elem_names[g] << " dimension "
                   << "(" << (gridded_tensor[6].*f[g])() << ") "
                   << "of the data tensor must be 1"
                   << endl << endl;
            }
        }
    }

}

