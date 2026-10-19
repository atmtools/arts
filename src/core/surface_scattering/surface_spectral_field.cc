
#include "surface_spectral_field.h"

#include <functional_field_interp.h>
#include <xml_io_base.h>

#include <algorithm>

namespace surface_scattering {
Data::Data()                            = default;
Data::Data(const Data &)                = default;
Data::Data(Data &&) noexcept            = default;
Data &Data::operator=(const Data &)     = default;
Data &Data::operator=(Data &&) noexcept = default;

void Data::adjust_interpolation_extrapolation() {
  using enum InterpolationExtrapolation;
  const auto [frq, lat, lon] = field_interp::single_point_dims(data);
  if (frq) frq_low = frq_upp = Nearest;
  if (lat) lat_low = lat_upp = Nearest;
  if (lon) lon_low = lon_upp = Nearest;
}

// Allow copy implicitly from all types
Data::Data(Numeric x) : data(x) { adjust_interpolation_extrapolation(); }

Data &Data::operator=(Numeric x) {
  data = x;
  adjust_interpolation_extrapolation();
  return *this;
}

Data::Data(SpectralGeodeticField3 x) : data(std::move(x)) { adjust_interpolation_extrapolation(); }

Data &Data::operator=(SpectralGeodeticField3 x) {
  data = std::move(x);
  adjust_interpolation_extrapolation();
  return *this;
}

Data::Data(FunctionalData x) : data(std::move(x)) { adjust_interpolation_extrapolation(); }

Data &Data::operator=(FunctionalData x) {
  data = std::move(x);
  adjust_interpolation_extrapolation();
  return *this;
}
String Data::data_type() const {
  if (std::holds_alternative<SpectralGeodeticField3>(data)) return "SpectralGeodeticField3";
  if (std::holds_alternative<Numeric>(data)) return "Numeric";
  if (std::holds_alternative<FunctionalData>(data)) return "FunctionalData";

  ARTS_USER_ERROR("Cannot understand data type; is this a new type")
}

ConstVectorView Data::flat_view() const { return field_interp::flat_view(data); }

VectorView Data::flat_view() { return field_interp::flat_view(data); }

//! Flat weights for the positions on the surface
std::vector<std::pair<Index, Numeric>> Data::flat_weights(const Numeric &frq,
                                                          const Numeric &lat,
                                                          const Numeric &lon) const {
  return field_interp::flat_weight(data, {frq, lat, lon});
}

Numeric Data::at(const Numeric frq, const Numeric lat, const Numeric lon) const {
  return field_interp::at(data, {{{frq_low, frq_upp}, {lat_low, lat_upp}, {lon_low, lon_upp}}}, {frq, lat, lon});
}

Vector Data::at(const Vector &frq, const Numeric lat, const Numeric lon) const {
  Vector out(frq.size());
  std::transform(frq.begin(), frq.end(), out.begin(), [&](Numeric f) { return at(f, lat, lon); });
  return out;
}

bool Data::ok() const { return field_interp::ok(data); }
}  //namespace surface_scattering

void xml_io_stream<SurfaceSpectralData>::write(std::ostream              &os,
                                               const SurfaceSpectralData &x,
                                               bofstream                 *pbofs,
                                               std::string_view           name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.data, pbofs, "Data"sv);
  xml_write_to_stream(os, x.frq_upp, pbofs, "frq_upp"sv);
  xml_write_to_stream(os, x.frq_low, pbofs, "frq_low"sv);
  xml_write_to_stream(os, x.lat_upp, pbofs, "lat_upp"sv);
  xml_write_to_stream(os, x.lat_low, pbofs, "lat_low"sv);
  xml_write_to_stream(os, x.lon_upp, pbofs, "lon_upp"sv);
  xml_write_to_stream(os, x.lon_low, pbofs, "lon_low"sv);

  tag.write_to_end_stream(os);
}

void xml_io_stream<SurfaceSpectralData>::read(std::istream &is, SurfaceSpectralData &x, bifstream *pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.data, pbifs);
  xml_read_from_stream(is, x.frq_upp, pbifs);
  xml_read_from_stream(is, x.frq_low, pbifs);
  xml_read_from_stream(is, x.lat_upp, pbifs);
  xml_read_from_stream(is, x.lat_low, pbifs);
  xml_read_from_stream(is, x.lon_upp, pbifs);
  xml_read_from_stream(is, x.lon_low, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
