#ifndef predefined_predef_data_h
#define predefined_predef_data_h

#include <debug.h>
#include <enums.h>

#include <exception>
#include <map>
#include <memory>
#include <ostream>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace Absorption::PredefinedModel {
  ENUMCLASS(DataKey, char, 
  water_mt_ckd_4d0
  )
namespace MT_CKD400 {
struct WaterData {
  static constexpr DataKey key = DataKey::water_mt_ckd_4d0;
  double ref_press;
  double ref_temp;
  double ref_h2o_vmr;
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
}  // namespace MT_CKD400

struct Model {
  using DataHolder = std::variant<std::monostate, MT_CKD400::WaterData>;
  using DataMap = std::map<DataKey, DataHolder>;

  DataMap data;
  
  Model() = default;
  Model(const Model&);
  Model& operator=(const Model&);

  template <typename T>
  [[nodiscard]] const T& get() const try {
    ARTS_USER_ERROR_IF(not good_enum(T::key), "Bad key")
    ARTS_USER_ERROR_IF(data.find(T::key) == data.end(), "No data")
    return std::get<T>(data.at(T::key));
  } catch (std::exception& e) {
    ARTS_USER_ERROR(
        "Cannot find data for continua with error: ",
        e.what())
  }

  template <typename T>
  void set(T x) {
    ARTS_USER_ERROR_IF(not good_enum(T::key), "Bad key")
    data[T::key] = std::move(x);
  }

  void set_all(const Model&);

  [[nodiscard]] std::vector<std::size_t> data_size(DataKey) const;
  [[nodiscard]] std::size_t size() const {return data.size();}
  void resize(const std::vector<std::size_t>&, DataKey);
  [[nodiscard]] std::vector<DataKey> keys() const;
  void output_data_to_stream(std::ostream&, DataKey) const;
  void set_data_from_stream(std::istream&, DataKey);

  friend std::ostream& operator<<(std::ostream&, const Model&);
};
}  // namespace Absorption::PredefinedModel

using PredefinedModelData = Absorption::PredefinedModel::Model;
using PredefinedModelDataKey = Absorption::PredefinedModel::DataKey;

#endif  // predefined_predef_data_h
