#ifndef predefined_predef_data_h
#define predefined_predef_data_h

#include <debug.h>

#include <exception>
#include <map>
#include <memory>
#include <ostream>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace Absorption::PredefinedModel {
namespace Hitran::MTCKD {
constexpr std::string_view water_key = "MTCKDWATER";
struct WaterData {
  std::vector<double> self_absco_ref;
  std::vector<double> for_absco_ref;
  std::vector<double> wavenumbers;
  std::vector<double> self_texp;
  friend std::ostream& operator<<(std::ostream&, const WaterData&);
  friend std::istream& operator>>(std::istream&, WaterData&);
  void resize(const std::vector<std::size_t>&);
};
}  // namespace Hitran::MTCKD

class Model {
public:
  using DataHolder = std::variant<std::unique_ptr<Hitran::MTCKD::WaterData>>;

private:
  using DataMap = std::map<std::string_view, DataHolder>;
  DataMap data;

public:
  Model() = default;
  Model(const Model&);
  Model& operator=(const Model&);

  template <typename T> [[nodiscard]] const T& get() const;
  void set(Hitran::MTCKD::WaterData x);

  [[nodiscard]] std::vector<std::string> keys() const;
  void output_data_to_stream(std::ostream&, std::string_view) const;
  void set_data_from_stream(std::istream&, std::string_view);

  friend std::ostream& operator<<(std::ostream&, const Model&);
};
}  // namespace Absorption::PredefinedModel

using PredefinedModelData = Absorption::PredefinedModel::Model;

#endif  // predefined_predef_data_h
