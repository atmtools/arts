#include "matpack.h"
#include "gridded_fields.h"

using std::cout;
using std::endl;


int main (void)
{
  // Creating two gridded fields
  //////////////////////////////////////////////////////////////////////////////
  GriddedField1 gfone("I'm a GriddedField1");
  GriddedField2 gftwo;
  
  gftwo.set_name ("I'm a GriddedField2");
  
  
  // Initializing the grids
  //////////////////////////////////////////////////////////////////////////////
  Vector gfonegrid(1,5,1);        // gfonegrid = [1,2,3,4,5]
  gfone.set_grid(0, gfonegrid);   // Set grid for the vector elements.
  
  Vector gftwogrid0(1,5,1);       // gftwogrid0 = [1,2,3,4,5]
  ArrayOfString gftwogrid1{"Chan1", "Chan2", "Chan3"};
  
  gftwo.set_grid(0, gftwogrid0);  // Set grid for the matrix rows.
  gftwo.set_grid(1, gftwogrid1);  // Set grid for the matrix columns.

  gfone.set_grid_name (0, "Pressure");
  
  gftwo.set_grid_name (0, "Pressure");
  gftwo.set_grid_name (1, "Channel");

  // Initializing the data
  //////////////////////////////////////////////////////////////////////////////
  Vector avector(1,4,0.5);    // avector = [1,1.5,2,2.5]
  
  gfone.data = avector;
  
  cout << gfone;

  Matrix amatrix(5,3,4.);     // amatrix = [[4,4,4],[4,4,4],...]
  
  gftwo.data = amatrix;
  
  // Consistency check
  //////////////////////////////////////////////////////////////////////////////
  
  if (!gfone.checksize())
    cout << gfone.get_name() << ": Sizes of grid and data don't match" << endl;
  
  // This should fail!
  if (!gftwo.checksize())
    cout << gftwo.get_name() << ": Sizes of grids and data don't match" << endl;
  
  // Output
  //////////////////////////////////////////////////////////////////////////////
  
  cout << "GriddedField1: " << gfone << endl;
  cout << "GriddedField2: " << gftwo << endl;
  
  return 0;
}

