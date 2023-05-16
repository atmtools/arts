/*!
  \file   m_tessem.cc

  \brief  This file contains functions that are adapted from TESSEM
  code which is used to calculate surface emissivity.
*/

#include "file.h"
#include "matpack_data.h"
#include "mystring.h"
#include "tessem.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void TessemNNReadAscii(TessemNN& net,
                       const String& net_file,
                       const Verbosity&) {
  ifstream net_is;

  open_input_file(net_is, net_file);
  tessem_read_ascii(net_is, net);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void TestTessem(Vector& outvalues,
                const TessemNN& net,
                const Vector& invalues,
                const Verbosity& verbosity) {
  CREATE_OUT1;
  outvalues.resize(net.nb_outputs);
  tessem_prop_nn(outvalues, net, invalues);
  out1 << "Input values    : " << invalues << "\n";
  out1 << "Output values   : " << outvalues << "\n";
}
