#include "matpack.h"
#include "make_array.h"
#include "gridded_fields.h"

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
  gfone.set_grid(0, gfonegrid);   // Set the grid for the first dimension.
  
  MakeArray<String> gftwogrid0("Channel 1", "Channel2", "Channel3");
  Vector gftwogrid1(1,5,1);       // gftwogrid1 = [1,2,3,4,5]
  
  gftwo.set_grid(0, gftwogrid0);  // Set grid for the first dimension
  gftwo.set_grid(1, gftwogrid1);  // Set grid for the second dimension

  gfone.set_grid_name (0, "Pressure");
  
  gftwo.set_grid_name (0, "Instrument channel");
  gftwo.set_grid_name (1, "Pressure");

  // Initializing the data
  //////////////////////////////////////////////////////////////////////////////
  Vector avector(1,5,0.5);    // avector = [1,1.5,2,2.5,3]
  
  (Vector&)gfone = avector;
  
  Matrix amatrix(2,3,4.);     // amatrix = [[4,4,4],[4,4,4]]
  
  (Matrix&)gftwo = amatrix;
  
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

