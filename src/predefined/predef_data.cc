#include "predef_data.h"

#include <memory>
#include <ostream>
#include <type_traits>

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

const Hitran::MTCKD::WaterData& Model::get_hitran_mtckd_water() const try {
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
}  // namespace Hitran::MTCKD

std::ostream& operator<<(std::ostream& os, const Model& data) {
  for (auto& x : data.data) {
    os << x.first << '\n';
    std::visit([&](auto& v) { os << (*v); }, x.second);
    os << '\n';
  }
  return os;
}
}  // namespace Absorption::PredefinedModel