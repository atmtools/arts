#include <memory>
#include "xml_io_scattering.h"

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::sht::SHT &sht,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  ArtsXMLTag open_tag, close_tag;
  open_tag.set_name("SHT");
  open_tag.add_attribute("l_max", sht.get_m_max());
  open_tag.add_attribute("m_max", sht.get_l_max());
  open_tag.add_attribute("n_aa", sht.get_n_azimuth_angles());
  open_tag.add_attribute("n_lat", sht.get_n_zenith_angles());
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

  Index l_max, m_max, n_za, n_lat;
  tag.get_attribute_value("l_max", l_max);
  tag.get_attribute_value("m_max", m_max);
  tag.get_attribute_value("n_za", n_za);
  tag.get_attribute_value("n_lat", n_lat);

  sht = scattering::sht::SHT(l_max, m_max, n_za, n_lat);
}



void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::IrregularZenithAngleGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &name) {
  ArtsXMLTag open_tag, close_tag;
  open_tag.set_name("IrregularZenithAngleGrid");
  open_tag.write_to_stream(os_xml);
  xml_write_to_stream(os_xml, static_cast<const Vector&>(grid), pbofs, name);
  close_tag.set_name("/IrregularZenithAngleGrid");
  close_tag.write_to_stream(os_xml);
  os_xml << '\n';
}

void xml_read_from_stream(std::istream &is_xml,
                          scattering::IrregularZenithAngleGrid& grid,
                          bifstream *pbifs [[maybe_unused]]) {
  ArtsXMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("IrregularZenithAngleGrid");
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

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::ZenithAngleGrid &grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &name) {
  ArtsXMLTag open_tag, close_tag;
  std::visit(
             [&](const auto& grd) {xml_write_to_stream(os_xml, grd, pbofs, name); },
             grid
             );
}

void xml_read_from_stream(std::istream &is_xml,
                          scattering::ZenithAngleGrid& za_grid,
                          bifstream *pbifs [[maybe_unused]]) {
  try {
    scattering::IrregularZenithAngleGrid il_grid{};
    xml_read_from_stream(is_xml, il_grid, pbifs);
    za_grid = il_grid;
    return;
  } catch (const std::runtime_error& e) {}

  try {
    scattering::GaussLegendreGrid gl_grid{};
    xml_read_from_stream(is_xml, gl_grid, pbifs);
    za_grid = gl_grid;
    return;
  } catch (const std::runtime_error& e) {}

  try {
    scattering::DoubleGaussGrid dg_grid{};
    xml_read_from_stream(is_xml, dg_grid, pbifs);
    za_grid = dg_grid;
    return;
  } catch (const std::runtime_error& e) {}

  try {
    scattering::LobattoGrid l_grid{};
    xml_read_from_stream(is_xml, l_grid, pbifs);
    za_grid = l_grid;
    return;
  } catch (const std::runtime_error& e) {}

  try {
    scattering::FejerGrid f_grid{};
    xml_read_from_stream(is_xml, f_grid, pbifs);
    za_grid = f_grid;
    return;
  } catch (const std::runtime_error& e) {}

  xml_parse_error("Encountered unknown zenith angle grid.");
}




template <std::floating_point Scalar, scattering::Format fmt, scattering::Representation repr>
void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::PhaseMatrixData<Scalar, fmt, repr>& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &) {
  //if constexpr (std::is_same_v<format, scattering::Format::ARO>) {
  //  ArtsXMLTag open_tag, close_tag;
  //  open_tag.set_name("PhaseMatrixDataAROGridded");
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

template <std::floating_point Scalar, scattering::Format fmt, scattering::Representation repr>
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
