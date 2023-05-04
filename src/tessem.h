/*!
  \file   tessem.h

  \brief  This file contains functions that are adapted from TESSEM
  code which is used to calculate surface emissivity.
*/

#ifndef tessem_h
#define tessem_h

#include <fstream>
#include "matpack_data.h"

struct TessemNN {
  Index nb_inputs;
  Index nb_outputs;
  Index nb_cache;
  Vector b1;
  Vector b2;
  Matrix w1;
  Matrix w2;
  Vector x_min;
  Vector x_max;
  Vector y_min;
  Vector y_max;

  friend std::ostream& operator<<(std::ostream& os, const TessemNN&) {return os;}
};

void tessem_read_ascii(std::ifstream& is, TessemNN& net);

void tessem_prop_nn(VectorView ny, const TessemNN& net, ConstVectorView nx);

#endif /* tessem_h */
