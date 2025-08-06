#include "xml_atm.h"

void xml_io_stream<AtmPoint>::read(std::istream& is,
                                   AtmPoint& v,
                                   bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, v.specs, pbifs);
  xml_read_from_stream(is, v.isots, pbifs);
  xml_read_from_stream(is, v.nlte, pbifs);
  xml_read_from_stream(is, v.ssprops, pbifs);
  xml_read_from_stream(is, v.pressure, pbifs);
  xml_read_from_stream(is, v.temperature, pbifs);
  xml_read_from_stream(is, v.wind, pbifs);
  xml_read_from_stream(is, v.mag, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<AtmPoint>::write(std::ostream& os,
                                    const AtmPoint& v,
                                    bofstream* pbofs,
                                    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, v.specs, pbofs, "Species Data"sv);
  xml_write_to_stream(os, v.isots, pbofs, "Isotopologue Data"sv);
  xml_write_to_stream(os, v.nlte, pbofs, "NLTE Data"sv);
  xml_write_to_stream(os, v.ssprops, pbofs, "Scattering Data"sv);
  xml_write_to_stream(os, v.pressure, pbofs, "pressure"sv);
  xml_write_to_stream(os, v.temperature, pbofs, "temperature"sv);
  xml_write_to_stream(os, v.wind, pbofs, "wind field"sv);
  xml_write_to_stream(os, v.mag, pbofs, "magnetic field"sv);

  tag.write_to_end_stream(os);
}

void xml_io_stream<AtmData>::read(std::istream& is,
                                  AtmData& v,
                                  bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, v.data, pbifs);
  xml_read_from_stream(is, v.alt_upp, pbifs);
  xml_read_from_stream(is, v.alt_low, pbifs);
  xml_read_from_stream(is, v.lat_upp, pbifs);
  xml_read_from_stream(is, v.lat_low, pbifs);
  xml_read_from_stream(is, v.lon_upp, pbifs);
  xml_read_from_stream(is, v.lon_low, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<AtmData>::write(std::ostream& os,
                                   const AtmData& v,
                                   bofstream* pbofs,
                                   std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, v.data, pbofs, "Data"sv);
  xml_write_to_stream(os, v.alt_upp, pbofs, "alt_upp"sv);
  xml_write_to_stream(os, v.alt_low, pbofs, "alt_low"sv);
  xml_write_to_stream(os, v.lat_upp, pbofs, "lat_upp"sv);
  xml_write_to_stream(os, v.lat_low, pbofs, "lat_low"sv);
  xml_write_to_stream(os, v.lon_upp, pbofs, "lon_upp"sv);
  xml_write_to_stream(os, v.lon_low, pbofs, "lon_low"sv);

  tag.write_to_end_stream(os);
}

void xml_io_stream<AtmField>::read(std::istream& is,
                                   AtmField& v,
                                   bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, v.top_of_atmosphere, pbifs);

  xml_read_from_stream(is, v.other, pbifs);
  xml_read_from_stream(is, v.specs, pbifs);
  xml_read_from_stream(is, v.isots, pbifs);
  xml_read_from_stream(is, v.nlte, pbifs);
  xml_read_from_stream(is, v.ssprops, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<AtmField>::write(std::ostream& os,
                                    const AtmField& v,
                                    bofstream* pbofs,
                                    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, v.top_of_atmosphere, pbofs, "TOA"sv);

  xml_write_to_stream(os, v.other, pbofs, "Other"sv);
  xml_write_to_stream(os, v.specs, pbofs, "Species Data"sv);
  xml_write_to_stream(os, v.isots, pbofs, "Isotopologue Data"sv);
  xml_write_to_stream(os, v.nlte, pbofs, "NLTE Data"sv);
  xml_write_to_stream(os, v.ssprops, pbofs, "Scattering Data"sv);

  tag.write_to_end_stream(os);
}
