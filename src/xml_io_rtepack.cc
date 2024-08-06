#include <rtepack.h>

#include "xml_io.h"
#include "xml_io_array_macro.h"

//=== Propmat ================================================================

//! Reads Propmat from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      Propmat return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          Propmat &pm,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("Propmat");

  if (pbifs)
    pbifs->readDoubleArray(pm.data.data(), 7);
  else
    is_xml >> pm.A() >> pm.B() >> pm.C() >> pm.D() >> pm.U() >> pm.V() >>
        pm.W();

  tag.read_from_stream(is_xml);
  tag.check_name("/Propmat");
}

//! Writes Propmat to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      Propmat
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const Propmat &pm,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Propmat");
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  if (pbofs)
    *pbofs << pm.A() << pm.B() << pm.C() << pm.D() << pm.U() << pm.V()
           << pm.W();
  else
    os_xml << ' ' << pm.A() << ' ' << pm.B() << ' ' << pm.C() << ' ' << pm.D()
           << ' ' << pm.U() << ' ' << pm.V() << ' ' << pm.W() << ' ';

  close_tag.set_name("/Propmat");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== Stokvec ================================================================

//! Reads Stokvec from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      Stokvec return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          Stokvec &pm,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("Stokvec");

  if (pbifs)
    pbifs->readDoubleArray(pm.data.data(), 4);
  else
    is_xml >> pm.I() >> pm.Q() >> pm.U() >> pm.V();

  tag.read_from_stream(is_xml);
  tag.check_name("/Stokvec");
}

//! Writes Stokvec to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      Stokvec
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const Stokvec &pm,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Stokvec");
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  if (pbofs)
    *pbofs << pm.I() << pm.Q() << pm.U() << pm.V();
  else
    os_xml << ' ' << pm.I() << ' ' << pm.Q() << ' ' << pm.U() << ' ' << pm.V()
           << ' ';

  close_tag.set_name("/Stokvec");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== Muelmat ================================================================

//! Reads Muelmat from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      Muelmat return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          Muelmat &pm,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("Muelmat");

  if (pbifs)
    pbifs->readDoubleArray(pm.data.data(), 16);
  else
    is_xml >> pm.data[0] >> pm.data[1] >> pm.data[2] >> pm.data[3] >>
        pm.data[4] >> pm.data[5] >> pm.data[6] >> pm.data[7] >> pm.data[8] >>
        pm.data[9] >> pm.data[10] >> pm.data[11] >> pm.data[12] >>
        pm.data[13] >> pm.data[14] >> pm.data[15];

  tag.read_from_stream(is_xml);
  tag.check_name("/Muelmat");
}

//! Writes Muelmat to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      Muelmat
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const Muelmat &pm,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("Muelmat");
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  if (pbofs)
    *pbofs << pm.data[0] << pm.data[1] << pm.data[2] << pm.data[3] << pm.data[4]
           << pm.data[5] << pm.data[6] << pm.data[7] << pm.data[8] << pm.data[9]
           << pm.data[10] << pm.data[11] << pm.data[12] << pm.data[13]
           << pm.data[14] << pm.data[15];
  else
    os_xml << ' ' << pm.data[0] << ' ' << pm.data[1] << ' ' << pm.data[2] << ' '
           << pm.data[3] << '\n'
           << ' ' << pm.data[4] << ' ' << pm.data[5] << ' ' << pm.data[6] << ' '
           << pm.data[7] << '\n'
           << ' ' << pm.data[8] << ' ' << pm.data[9] << ' ' << pm.data[10]
           << ' ' << pm.data[11] << '\n'
           << ' ' << pm.data[12] << ' ' << pm.data[13] << ' ' << pm.data[14]
           << ' ' << pm.data[15] << ' ';

  close_tag.set_name("/Muelmat");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== PropmatVector ================================================================

//! Reads PropmatVector from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      PropmatVector return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          PropmatVector &pmv,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("PropmatVector");

  Index n;
  tag.get_attribute_value("nelem", n);

  pmv.resize(n);
  for (auto &pm : pmv) {
    if (pbifs)
      pbifs->readDoubleArray(pm.data.data(), 7);
    else
      is_xml >> pm.A() >> pm.B() >> pm.C() >> pm.D() >> pm.U() >> pm.V() >>
          pm.W();
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/PropmatVector");
}

//! Writes PropmatVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      PropmatVector
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const PropmatVector &pmv,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("PropmatVector");
  open_tag.add_attribute("nelem", pmv.nelem());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto &pm : pmv) {
    if (pbofs)
      *pbofs << pm.A() << pm.B() << pm.C() << pm.D() << pm.U() << pm.V()
             << pm.W();
    else
      os_xml << pm.A() << ' ' << pm.B() << ' ' << pm.C() << ' ' << pm.D() << ' '
             << pm.U() << ' ' << pm.V() << ' ' << pm.W() << '\n';
  }

  close_tag.set_name("/PropmatVector");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== StokvecVector ================================================================

//! Reads StokvecVector from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      StokvecVector return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          StokvecVector &pmv,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("StokvecVector");

  Index n;
  tag.get_attribute_value("nelem", n);

  pmv.resize(n);
  for (auto &pm : pmv) {
    if (pbifs)
      pbifs->readDoubleArray(pm.data.data(), 4);
    else
      is_xml >> pm.I() >> pm.Q() >> pm.U() >> pm.V();
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/StokvecVector");
}

//! Writes StokvecVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      StokvecVector
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const StokvecVector &pmv,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("StokvecVector");
  open_tag.add_attribute("nelem", pmv.nelem());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto &pm : pmv) {
    if (pbofs)
      *pbofs << pm.I() << pm.Q() << pm.U() << pm.V();
    else
      os_xml << ' ' << pm.I() << ' ' << pm.Q() << ' ' << pm.U() << ' ' << pm.V()
             << '\n';
  }

  close_tag.set_name("/StokvecVector");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== MuelmatVector ================================================================

//! Reads MuelmatVector from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      MuelmatVector return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          MuelmatVector &pmv,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("MuelmatVector");

  Index n;
  tag.get_attribute_value("nelem", n);

  pmv.resize(n);
  for (auto &pm : pmv) {
    if (pbifs)
      pbifs->readDoubleArray(pm.data.data(), 16);
    else
      is_xml >> pm.data[0] >> pm.data[1] >> pm.data[2] >> pm.data[3] >>
          pm.data[4] >> pm.data[5] >> pm.data[6] >> pm.data[7] >> pm.data[8] >>
          pm.data[9] >> pm.data[10] >> pm.data[11] >> pm.data[12] >>
          pm.data[13] >> pm.data[14] >> pm.data[15];
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/MuelmatVector");
}

//! Writes MuelmatVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      MuelmatVector
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const MuelmatVector &pmv,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("MuelmatVector");
  open_tag.add_attribute("nelem", pmv.nelem());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto &pm : pmv) {
    if (pbofs)
      *pbofs << pm.data[0] << pm.data[1] << pm.data[2] << pm.data[3]
             << pm.data[4] << pm.data[5] << pm.data[6] << pm.data[7]
             << pm.data[8] << pm.data[9] << pm.data[10] << pm.data[11]
             << pm.data[12] << pm.data[13] << pm.data[14] << pm.data[15];
    else
      os_xml << ' ' << pm.data[0] << ' ' << pm.data[1] << ' ' << pm.data[2]
             << ' ' << pm.data[3] << '\n'
             << ' ' << pm.data[4] << ' ' << pm.data[5] << ' ' << pm.data[6]
             << ' ' << pm.data[7] << '\n'
             << ' ' << pm.data[8] << ' ' << pm.data[9] << ' ' << pm.data[10]
             << ' ' << pm.data[11] << '\n'
             << ' ' << pm.data[12] << ' ' << pm.data[13] << ' ' << pm.data[14]
             << ' ' << pm.data[15] << '\n';
  }

  close_tag.set_name("/MuelmatVector");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== PropmatMatrix ================================================================

//! Reads PropmatMatrix from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      PropmatVector return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          PropmatMatrix &pmm,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("PropmatMatrix");

  Index nr, nc;
  tag.get_attribute_value("nrows", nr);
  tag.get_attribute_value("ncols", nc);

  pmm.resize(nr, nc);
  for (auto &&pmv : pmm) {
    for (auto &pm : pmv) {
      if (pbifs)
        pbifs->readDoubleArray(pm.data.data(), 7);
      else
        is_xml >> pm.A() >> pm.B() >> pm.C() >> pm.D() >> pm.U() >> pm.V() >>
            pm.W();
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/PropmatMatrix");
}

//! Writes PropmatMatrix to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      PropmatMatrix
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const PropmatMatrix &pmm,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("PropmatMatrix");
  open_tag.add_attribute("nrows", pmm.nrows());
  open_tag.add_attribute("ncols", pmm.ncols());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto &pmv : pmm) {
    for (const auto &pm : pmv) {
      if (pbofs)
        *pbofs << pm.A() << pm.B() << pm.C() << pm.D() << pm.U() << pm.V()
               << pm.W();
      else
        os_xml << pm.A() << ' ' << pm.B() << ' ' << pm.C() << ' ' << pm.D()
               << ' ' << pm.U() << ' ' << pm.V() << ' ' << pm.W() << '\n';
    }
  }

  close_tag.set_name("/PropmatMatrix");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== StokvecMatrix ================================================================

//! Reads StokvecMatrix from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      StokvecMatrix return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          StokvecMatrix &pmm,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("StokvecMatrix");

  Index nr, nc;
  tag.get_attribute_value("nrows", nr);
  tag.get_attribute_value("ncols", nc);

  pmm.resize(nr, nc);
  for (auto &&pmv : pmm) {
    for (auto &pm : pmv) {
      if (pbifs)
        pbifs->readDoubleArray(pm.data.data(), 4);
      else
        is_xml >> pm.I() >> pm.Q() >> pm.U() >> pm.V();
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/StokvecMatrix");
}

//! Writes StokvecMatrix to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      StokvecMatrix
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const StokvecMatrix &pmm,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("StokvecMatrix");
  open_tag.add_attribute("nrows", pmm.nrows());
  open_tag.add_attribute("ncols", pmm.ncols());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto &pmv : pmm) {
    for (const auto &pm : pmv) {
      if (pbofs)
        *pbofs << pm.I() << pm.Q() << pm.U() << pm.V();
      else
        os_xml << ' ' << pm.I() << ' ' << pm.Q() << ' ' << pm.U() << ' '
               << pm.V() << '\n';
    }
  }

  close_tag.set_name("/StokvecMatrix");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== StokvecTensor3 ================================================================

//! Reads StokvecTensor3 from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      StokvecTensor3 return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          StokvecTensor3 &pmt3,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("StokvecTensor3");

  Index np, nr, nc;
  tag.get_attribute_value("npages", np);
  tag.get_attribute_value("nrows", nr);
  tag.get_attribute_value("ncols", nc);

  pmt3.resize(np, nr, nc);
  for (auto &&pmm : pmt3) {
    for (auto &&pmv : pmm) {
      for (auto &pm : pmv) {
        if (pbifs)
          pbifs->readDoubleArray(pm.data.data(), 4);
        else
          is_xml >> pm.I() >> pm.Q() >> pm.U() >> pm.V();
      }
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/StokvecTensor3");
}

//! Writes StokvecTensor3 to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      StokvecTensor3
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const StokvecTensor3 &pmt3,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("StokvecTensor3");
  open_tag.add_attribute("npages", pmt3.npages());
  open_tag.add_attribute("nrows", pmt3.nrows());
  open_tag.add_attribute("ncols", pmt3.ncols());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  for (const auto &pmm : pmt3) {
    for (const auto &pmv : pmm) {
      for (const auto &pm : pmv) {
        if (pbofs)
          *pbofs << pm.I() << pm.Q() << pm.U() << pm.V();
        else
          os_xml << ' ' << pm.I() << ' ' << pm.Q() << ' ' << pm.U() << ' '
                 << pm.V() << '\n';
      }
    }
  }

  close_tag.set_name("/StokvecTensor3");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== StokvecTensor4 ================================================================

//! Reads StokvecTensor4 from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      StokvecTensor4 return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          StokvecTensor4 &pmt4,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("StokvecTensor4");

  Index nb, np, nr, nc;
  tag.get_attribute_value("nbooks", nb);
  tag.get_attribute_value("npages", np);
  tag.get_attribute_value("nrows", nr);
  tag.get_attribute_value("ncols", nc);

  pmt4.resize(nb, np, nr, nc);
  for (auto &&pmt3 : pmt4) {
    for (auto &&pmm : pmt3) {
      for (auto &&pmv : pmm) {
        for (auto &pm : pmv) {
          if (pbifs)
            pbifs->readDoubleArray(pm.data.data(), 4);
          else
            is_xml >> pm.I() >> pm.Q() >> pm.U() >> pm.V();
        }
      }
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/StokvecTensor4");
}

//! Writes StokvecTensor4 to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      StokvecTensor4
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const StokvecTensor4 &pmt4,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("StokvecTensor4");
  open_tag.add_attribute("nbooks", pmt4.nbooks());
  open_tag.add_attribute("npages", pmt4.npages());
  open_tag.add_attribute("nrows", pmt4.nrows());
  open_tag.add_attribute("ncols", pmt4.ncols());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  for (const auto &pmt3 : pmt4) {
    for (const auto &pmm : pmt3) {
      for (const auto &pmv : pmm) {
        for (const auto &pm : pmv) {
          if (pbofs)
            *pbofs << pm.I() << pm.Q() << pm.U() << pm.V();
          else
            os_xml << ' ' << pm.I() << ' ' << pm.Q() << ' ' << pm.U() << ' '
                   << pm.V() << '\n';
        }
      }
    }
  }

  close_tag.set_name("/StokvecTensor4");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== StokvecTensor5 ================================================================

//! Reads StokvecTensor5 from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      StokvecTensor5 return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          StokvecTensor5 &pmt5,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("StokvecTensor5");

  Index nb, np, nr, nc, ns;
  tag.get_attribute_value("nshelves", ns);
  tag.get_attribute_value("nbooks", nb);
  tag.get_attribute_value("npages", np);
  tag.get_attribute_value("nrows", nr);
  tag.get_attribute_value("ncols", nc);

  pmt5.resize(ns, nb, np, nr, nc);
  for (auto &&pmt4 : pmt5) {
    for (auto &&pmt3 : pmt4) {
      for (auto &&pmm : pmt3) {
        for (auto &&pmv : pmm) {
          for (auto &pm : pmv) {
            if (pbifs)
              pbifs->readDoubleArray(pm.data.data(), 4);
            else
              is_xml >> pm.I() >> pm.Q() >> pm.U() >> pm.V();
          }
        }
      }
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/StokvecTensor5");
}

//! Writes StokvecTensor5 to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      StokvecTensor5
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const StokvecTensor5 &pmt5,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("StokvecTensor5");
  open_tag.add_attribute("nshelves", pmt5.nshelves());
  open_tag.add_attribute("nbooks", pmt5.nbooks());
  open_tag.add_attribute("npages", pmt5.npages());
  open_tag.add_attribute("nrows", pmt5.nrows());
  open_tag.add_attribute("ncols", pmt5.ncols());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  for (const auto &pmt4 : pmt5) {
    for (const auto &pmt3 : pmt4) {
      for (const auto &pmm : pmt3) {
        for (const auto &pmv : pmm) {
          for (const auto &pm : pmv) {
            if (pbofs)
              *pbofs << pm.I() << pm.Q() << pm.U() << pm.V();
            else
              os_xml << ' ' << pm.I() << ' ' << pm.Q() << ' ' << pm.U() << ' '
                     << pm.V() << '\n';
          }
        }
      }
    }
  }

  close_tag.set_name("/StokvecTensor5");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== StokvecTensor6 ================================================================

//! Reads StokvecTensor6 from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      StokvecTensor6 return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          StokvecTensor6 &pmt6,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("StokvecTensor6");

  Index nb, np, nr, nc, ns, nv;
  tag.get_attribute_value("nvitrines", nv);
  tag.get_attribute_value("nshelves", ns);
  tag.get_attribute_value("nbooks", nb);
  tag.get_attribute_value("npages", np);
  tag.get_attribute_value("nrows", nr);
  tag.get_attribute_value("ncols", nc);

  pmt6.resize(nv, ns, nb, np, nr, nc);
  for (auto &&pmt5 : pmt6) {
    for (auto &&pmt4 : pmt5) {
      for (auto &&pmt3 : pmt4) {
        for (auto &&pmm : pmt3) {
          for (auto &&pmv : pmm) {
            for (auto &pm : pmv) {
              if (pbifs)
                pbifs->readDoubleArray(pm.data.data(), 4);
              else
                is_xml >> pm.I() >> pm.Q() >> pm.U() >> pm.V();
            }
          }
        }
      }
    }
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/StokvecTensor6");
}

//! Writes StokvecTensor6 to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      StokvecTensor6
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const StokvecTensor6 &pmt6,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("StokvecTensor6");
  open_tag.add_attribute("nvitrines", pmt6.nvitrines());
  open_tag.add_attribute("nshelves", pmt6.nshelves());
  open_tag.add_attribute("nbooks", pmt6.nbooks());
  open_tag.add_attribute("npages", pmt6.npages());
  open_tag.add_attribute("nrows", pmt6.nrows());
  open_tag.add_attribute("ncols", pmt6.ncols());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);
  for (const auto &pmt5 : pmt6) {
    for (const auto &pmt4 : pmt5) {
      for (const auto &pmt3 : pmt4) {
        for (const auto &pmm : pmt3) {
          for (const auto &pmv : pmm) {
            for (const auto &pm : pmv) {
              if (pbofs)
                *pbofs << pm.I() << pm.Q() << pm.U() << pm.V();
              else
                os_xml << ' ' << pm.I() << ' ' << pm.Q() << ' ' << pm.U() << ' '
                       << pm.V() << '\n';
            }
          }
        }
      }
    }
  }

  close_tag.set_name("/StokvecTensor6");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== MuelmatMatrix ================================================================

//! Reads MuelmatMatrix from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      MuelmatMatrix return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          MuelmatMatrix &pmm,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("MuelmatMatrix");

  Index nr, nc;
  tag.get_attribute_value("nrows", nr);
  tag.get_attribute_value("ncols", nc);

  pmm.resize(nr, nc);
  for (auto &&pmv : pmm) {
    for (auto &pm : pmv) {
      if (pbifs)
        pbifs->readDoubleArray(pm.data.data(), 16);
      else
        is_xml >> pm.data[0] >> pm.data[1] >> pm.data[2] >> pm.data[3] >>
            pm.data[4] >> pm.data[5] >> pm.data[6] >> pm.data[7] >>
            pm.data[8] >> pm.data[9] >> pm.data[10] >> pm.data[11] >>
            pm.data[12] >> pm.data[13] >> pm.data[14] >> pm.data[15];
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/MuelmatMatrix");
}

//! Writes MuelmatMatrix to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      MuelmatMatrix
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const MuelmatMatrix &pmm,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("MuelmatMatrix");
  open_tag.add_attribute("nrows", pmm.nrows());
  open_tag.add_attribute("ncols", pmm.ncols());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto &pmv : pmm) {
    for (const auto &pm : pmv) {
      if (pbofs)
        *pbofs << pm.data[0] << pm.data[1] << pm.data[2] << pm.data[3]
               << pm.data[4] << pm.data[5] << pm.data[6] << pm.data[7]
               << pm.data[8] << pm.data[9] << pm.data[10] << pm.data[11]
               << pm.data[12] << pm.data[13] << pm.data[14] << pm.data[15];
      else
        os_xml << ' ' << pm.data[0] << ' ' << pm.data[1] << ' ' << pm.data[2]
               << ' ' << pm.data[3] << '\n'
               << ' ' << pm.data[4] << ' ' << pm.data[5] << ' ' << pm.data[6]
               << ' ' << pm.data[7] << '\n'
               << ' ' << pm.data[8] << ' ' << pm.data[9] << ' ' << pm.data[10]
               << ' ' << pm.data[11] << '\n'
               << ' ' << pm.data[12] << ' ' << pm.data[13] << ' ' << pm.data[14]
               << ' ' << pm.data[15] << '\n';
    }
  }

  close_tag.set_name("/MuelmatMatrix");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

//=== MuelmatTensor3 ================================================================

//! Reads MuelmatTensor3 from XML input stream
/*!
 *  \param is_xml  XML Input stream
 *  \param pm      MuelmatTensor3 return value
 *  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(std::istream &is_xml,
                          MuelmatTensor3 &pmt3,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("MuelmatTensor3");

  Index nr, nc, np;
  tag.get_attribute_value("npages", np);
  tag.get_attribute_value("nrows", nr);
  tag.get_attribute_value("ncols", nc);

  pmt3.resize(np, nr, nc);
  for (auto &&pmm : pmt3) {
    for (auto &&pmv : pmm) {
      for (auto &pm : pmv) {
        if (pbifs)
          pbifs->readDoubleArray(pm.data.data(), 16);
        else
          is_xml >> pm.data[0] >> pm.data[1] >> pm.data[2] >> pm.data[3] >>
              pm.data[4] >> pm.data[5] >> pm.data[6] >> pm.data[7] >>
              pm.data[8] >> pm.data[9] >> pm.data[10] >> pm.data[11] >>
              pm.data[12] >> pm.data[13] >> pm.data[14] >> pm.data[15];
      }
    }
  }
  tag.read_from_stream(is_xml);
  tag.check_name("/MuelmatTensor3");
}

//! Writes MuelmatTensor3 to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param pm      MuelmatTensor3
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */
void xml_write_to_stream(std::ostream &os_xml,
                         const MuelmatTensor3 &pmt3,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name("MuelmatTensor3");
  open_tag.add_attribute("npages", pmt3.npages());
  open_tag.add_attribute("nrows", pmt3.nrows());
  open_tag.add_attribute("ncols", pmt3.ncols());
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  xml_set_stream_precision(os_xml);

  for (const auto &pmm : pmt3) {
    for (const auto &pmv : pmm) {
      for (const auto &pm : pmv) {
        if (pbofs)
          *pbofs << pm.data[0] << pm.data[1] << pm.data[2] << pm.data[3]
                 << pm.data[4] << pm.data[5] << pm.data[6] << pm.data[7]
                 << pm.data[8] << pm.data[9] << pm.data[10] << pm.data[11]
                 << pm.data[12] << pm.data[13] << pm.data[14] << pm.data[15];
        else
          os_xml << ' ' << pm.data[0] << ' ' << pm.data[1] << ' ' << pm.data[2]
                 << ' ' << pm.data[3] << '\n'
                 << ' ' << pm.data[4] << ' ' << pm.data[5] << ' ' << pm.data[6]
                 << ' ' << pm.data[7] << '\n'
                 << ' ' << pm.data[8] << ' ' << pm.data[9] << ' ' << pm.data[10]
                 << ' ' << pm.data[11] << '\n'
                 << ' ' << pm.data[12] << ' ' << pm.data[13] << ' '
                 << pm.data[14] << ' ' << pm.data[15] << '\n';
      }
    }
  }

  close_tag.set_name("/MuelmatTensor3");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

// arrays
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfPropmatVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfPropmatMatrix)

TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfStokvecVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfStokvecMatrix)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfStokvecTensor3)

TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfMuelmatVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfMuelmatMatrix)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfMuelmatTensor3)

TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfPropmatVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfPropmatMatrix)

TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfStokvecVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfStokvecMatrix)

TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfMuelmatVector)
TMPL_XML_READ_WRITE_STREAM_ARRAY(ArrayOfArrayOfMuelmatMatrix)
