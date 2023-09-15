/**
 * @file   ppath_struct.cc
 * @author Manfred Brath <manfred.brath@uni-hamburg.de>
 * @date   2022-09-028
 *
 * @brief  Propagation path structure functions.
 *
 * This file contains the definition of the Ppath structure function.
 */

#include "ppath_struct.h"

std::ostream& operator<<(std::ostream& os, const Ppath& x) {
  os << "dim: " << x.dim << "\n";
  os << "np: " << x.np << "\n";
  //os << "constant: " << x.constant << "\n";   // ZZZ
  os << "background: " << x.background << "\n";  // ZZZ
  os << "start_pos: " << x.start_pos << "\n";
  os << "start_los: " << x.start_los << "\n";
  os << "start_lstep: " << x.start_lstep << "\n";
  os << "pos: " << x.pos << "\n";
  os << "los: " << x.los << "\n";
  // os << "r: " << x.r << "\n";   // ZZZ
  os << "lstep: " << x.lstep << "\n";
  os << "end_pos: " << x.end_pos << "\n";
  os << "end_los: " << x.end_los << "\n";
  os << "end_lstep: " << x.end_lstep << "\n";
  os << "nreal: " << x.nreal << "\n";
  os << "ngroup: " << x.ngroup << "\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfPpath& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}
