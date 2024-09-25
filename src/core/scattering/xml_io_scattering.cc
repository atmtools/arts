#include "xml_io_scattering.h"

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::sht::SHT &sht,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag, close_tag;
  open_tag.set_name("SHT");
  open_tag.add_attribute("l_max", sht.get_m_max());
  open_tag.add_attribute("m_max", sht.get_l_max());
  open_tag.add_attribute("n_lon", sht.get_n_longitudes());
  open_tag.add_attribute("n_lat", sht.get_n_latitudes());
  open_tag.write_to_stream(os_xml);
  close_tag.set_name("/SHT");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

void xml_read_from_stream(std::istream &is_xml,
                          scattering::sht::SHT &sht,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("SHT");

  Index l_max, m_max, n_lon, n_lat;
  tag.get_attribute_value("l_max", l_max);
  tag.get_attribute_value("m_max", m_max);
  tag.get_attribute_value("n_lon", n_lon);
  tag.get_attribute_value("n_lat", n_lat);

  sht = scattering::sht::SHT(l_max, m_max, n_lon, n_lat);
}



void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::IrregularLatitudeGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &name) {
  ArtsXMLTag open_tag, close_tag;
  open_tag.set_name("IrregularLatitudeGrid");
  open_tag.write_to_stream(os_xml);
  xml_write_to_stream(os_xml, static_cast<const Vector&>(grid), pbofs, name);
  close_tag.set_name("/IrregularLatitudeGrid");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

void xml_read_from_stream(std::istream &is_xml,
                          scattering::IrregularLatitudeGrid& grid,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("IrregularLatitudeGrid");
  xml_read_from_stream(is_xml, static_cast<Vector&>(grid), pbifs);
}


void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::GaussLegendreGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag, close_tag;
  open_tag.set_name("GaussLegendreGrid");
  open_tag.add_attribute("n", grid.get_degree());
  open_tag.write_to_stream(os_xml);
  close_tag.set_name("/GaussLegendreGrid");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

void xml_read_from_stream(std::istream &is_xml,
                          scattering::GaussLegendreGrid& grid,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("GaussLegendreGrid");

  Index n;
  tag.get_attribute_value("n", n);
  grid = scattering::GaussLegendreGrid(n);
}

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::DoubleGaussGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag, close_tag;
  open_tag.set_name("DoubleGaussGrid");
  open_tag.add_attribute("n", grid.get_degree());
  open_tag.write_to_stream(os_xml);
  close_tag.set_name("/DoubleGaussGrid");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

void xml_read_from_stream(std::istream &is_xml,
                          scattering::DoubleGaussGrid& grid,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("DoubleGaussGrid");

  Index n;
  tag.get_attribute_value("n", n);
  grid = scattering::DoubleGaussGrid(n);
}

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::LobattoGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag, close_tag;
  open_tag.set_name("LobattoGrid");
  open_tag.add_attribute("n", grid.get_degree());
  open_tag.write_to_stream(os_xml);
  close_tag.set_name("/LobattoGrid");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

void xml_read_from_stream(std::istream &is_xml,
                          scattering::LobattoGrid& grid,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("LobattoGrid");

  Index n;
  tag.get_attribute_value("n", n);
  grid = scattering::LobattoGrid(n);
}

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::FejerGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag, close_tag;
  open_tag.set_name("FejerGrid");
  open_tag.add_attribute("n", grid.get_degree());
  open_tag.write_to_stream(os_xml);
  close_tag.set_name("/FejerGrid");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

void xml_read_from_stream(std::istream &is_xml,
                          scattering::FejerGrid& grid,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("FejerGrid");

  Index n;
  tag.get_attribute_value("n", n);
  grid = scattering::FejerGrid(n);
}


template <std::floating_point Scalar, scattering::Format fmt, scattering::Representation repr, Index stokes_dim>
void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::PhaseMatrixData<Scalar, fmt, repr, stokes_dim>& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  //if constexpr (std::is_same_v<format, scattering::Format::ARO>) {
  //  ArtsXMLTag open_tag, close_tag;
  //  open_tag.set_name("PhaseMatrixDataAROGridded");
  //  open_tag.add_attribute("stokes_dim", stokes_dim);
  //  open_tag.write_to_stream(os_xml);
  //  close_tag.set_name("/FejerGrid");
  //  close_tag.write_to_stream(os_xml);
  //os_xml << '\n';
  //  if constexpr (std::is_same_v<repr, scattering::Format::Gridded>) {

  //  } else {

  //  }
  //} else {
  //  if constexpr (std::is_same_v<repr, scattering::Format::Gridded>) {

  //  } else {

  //  }
  //
  //}
}
