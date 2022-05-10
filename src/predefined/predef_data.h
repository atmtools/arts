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

#include <enums.h>

namespace Absorption::PredefinedModel {
  ENUMCLASS(PredefinedModelDataKey, char, 
  HITRANMTCKDWATER
  )
namespace Hitran::MTCKD {
struct WaterData {
  static constexpr PredefinedModelDataKey key = PredefinedModelDataKey::HITRANMTCKDWATER;
  std::vector<double> self_absco_ref;
  std::vector<double> for_absco_ref;
  std::vector<double> wavenumbers;
  std::vector<double> self_texp;

  void resize(const std::vector<std::size_t>&);
  [[nodiscard]] std::vector<std::size_t> sizes() const {
    return {self_absco_ref.size()};
  };

  friend std::ostream& operator<<(std::ostream&, const WaterData&);
  friend std::istream& operator>>(std::istream&, WaterData&);
};
}  // namespace Hitran::MTCKD

class Model {
public:
  using DataHolder = std::variant<std::unique_ptr<Hitran::MTCKD::WaterData>>;

private:
  using DataMap = std::map<PredefinedModelDataKey, DataHolder>;
  DataMap data;

public:
  Model() = default;
  Model(const Model&);
  Model& operator=(const Model&);

  template <typename T>
  [[nodiscard]] const T& get() const try {
    ARTS_USER_ERROR_IF(not good_enum(T::key), "Bad key")
    return *std::get<std::unique_ptr<T>>(data.at(T::key));
  } catch (std::exception& e) {
    ARTS_USER_ERROR(
        "Cannot find data for Hitran MT CKD Water continua with error: ",
        e.what())
  }

  template <typename T>
  void set(T x) {
    ARTS_USER_ERROR_IF(not good_enum(T::key), "Bad key")
    data[T::key] = std::make_unique<T>(std::move(x));
  }

  void set_all(const Model&);

  [[nodiscard]] std::vector<std::size_t> data_size(PredefinedModelDataKey) const;
  [[nodiscard]] std::size_t size() const {return data.size();}
  void resize(const std::vector<std::size_t>&, PredefinedModelDataKey);
  [[nodiscard]] std::vector<PredefinedModelDataKey> keys() const;
  void output_data_to_stream(std::ostream&, PredefinedModelDataKey) const;
  void set_data_from_stream(std::istream&, PredefinedModelDataKey);

  friend std::ostream& operator<<(std::ostream&, const Model&);
};
}  // namespace Absorption::PredefinedModel

using PredefinedModelData = Absorption::PredefinedModel::Model;

#endif  // predefined_predef_data_h
