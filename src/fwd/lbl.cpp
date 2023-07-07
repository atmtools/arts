#include "lbl.h"

#include <numeric>

Complex lbl::full::at(Numeric f) const {
  return std::transform_reduce(models.begin(),
                               models.end(),
                               Complex{},
                               std::plus<>{},
                               [f](const auto& m) {
                                 return std::visit(
                                     [f](auto& mod) { return mod.at(f); }, m);
                               });
}

std::size_t lbl::full::size() const {
  return std::transform_reduce(models.begin(),
                               models.end(),
                               std::size_t{},
                               std::plus<>{},
                               [](const auto& m) {
                                 return std::visit(
                                     [](const auto& mod) { return mod.size(); },
                                     m);
                               });
}

ComplexVector lbl::full::at(const Vector& f) const {
  ComplexVector out(f.size());
  at(out, f);
  return out;
}

void lbl::full::at(ComplexVector& out, const Vector& fs) const {
  std::transform(fs.begin(), fs.end(), out.begin(), [this](const auto& f) {
    return at(f);
  });
}

static_assert(lbl::bandable<lbl::full>,
              "lbl::mtckd::single is not representative of a band of lines");
