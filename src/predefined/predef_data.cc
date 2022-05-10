#include "predef_data.h"

#include <memory>
#include <ostream>
#include <type_traits>

#include "debug.h"
#include "double_imanip.h"

namespace Absorption::PredefinedModel {
Model::Model(const Model& d) {
  for (auto& x : d.data) {
    std::visit(
        [&](auto&& v) {
          using MyType = std::remove_reference_t<decltype(*v)>;
          MyType y{*v};
          data[x.first] = std::make_unique<MyType>(std::move(y));
        },
        x.second);
  }
}

Model& Model::operator=(const Model& d) {
  for (auto& x : d.data) {
    std::visit(
        [&](auto&& v) {
          using MyType = std::remove_reference_t<decltype(*v)>;
          MyType y{*v};
          data[x.first] = std::make_unique<MyType>(std::move(y));
        },
        x.second);
  }
  return *this;
}

template <>
const Hitran::MTCKD::WaterData& Model::get() const try {
  return *std::get<std::unique_ptr<Hitran::MTCKD::WaterData>>(
      data.at(Hitran::MTCKD::water_key));
} catch (std::exception& e) {
  ARTS_USER_ERROR(
      "Cannot find data for Hitran MT CKD Water continua with error: ",
      e.what())
}

void Model::set(Hitran::MTCKD::WaterData x) {
  data[Hitran::MTCKD::water_key] =
      std::make_unique<Hitran::MTCKD::WaterData>(std::move(x));
}

namespace Hitran::MTCKD {
std::ostream& operator<<(std::ostream& os,
                         const Hitran::MTCKD::WaterData& data) {
  for (auto& x : data.for_absco_ref) os << x << ' ';
  os << '\n';
  for (auto& x : data.self_absco_ref) os << x << ' ';
  os << '\n';
  for (auto& x : data.wavenumbers) os << x << ' ';
  os << '\n';
  for (auto& x : data.self_texp) os << x << ' ';
  return os;
}

std::istream& operator>>(std::istream& is, WaterData& data) {
  for (auto& x : data.for_absco_ref) is >> double_imanip() >> x;
  for (auto& x : data.self_absco_ref) is >> double_imanip() >> x;
  for (auto& x : data.wavenumbers) is >> double_imanip() >> x;
  for (auto& x : data.self_texp) is >> double_imanip() >> x;
  return is;
}

void WaterData::resize(const std::vector<std::size_t>& inds) {
  ARTS_USER_ERROR_IF(inds.size() not_eq 1, "Expects only one size")
  self_absco_ref.resize(inds.front());
  for_absco_ref.resize(inds.front());
  wavenumbers.resize(inds.front());
  self_texp.resize(inds.front());
}
}  // namespace Hitran::MTCKD

std::vector<std::string> Model::keys() const {
  std::vector<std::string> out;
  out.reserve(data.size());
  for (auto& x : data) out.emplace_back(x.first);
  return out;
}

void Model::output_data_to_stream(std::ostream& os, std::string_view key) const {
  std::visit([&](auto& v) { os << (*v); }, data.at(key));
}

void Model::set_data_from_stream(std::istream& is, std::string_view key) {
  std::visit([&](auto& v) { is >> (*v); }, data[key]);
}

std::ostream& operator<<(std::ostream& os, const Model& data) {
  for (auto& x : data.data) {
    os << x.first << '\n';
    data.output_data_to_stream(os, x.first);
    os << '\n';
  }
  return os;
}
}  // namespace Absorption::PredefinedModel