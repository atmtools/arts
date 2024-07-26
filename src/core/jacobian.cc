#include "jacobian.h"

#include <algorithm>

#include "compare.h"

namespace Jacobian {
std::ostream& operator<<(std::ostream& os, const AtmTarget& target) {
  return os << "Atmosphere key value: " << target.type << ", starting at "
            << target.x_start << " of size " << target.x_size;
}

void AtmTarget::update(AtmField& atm, const Vector& x) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_model(atm, type, x.slice(x_start, x_size));
}

void AtmTarget::update(Vector& x, const AtmField& atm) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_state(x.slice(x_start, x_size), atm, type);
}

std::ostream& operator<<(std::ostream& os, const SurfaceTarget& target) {
  return os << "Surface key value: " << target.type << ", starting at "
            << target.x_start << " of size " << target.x_size;
}

void SurfaceTarget::update(SurfaceField& surf, const Vector& x) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_model(surf, type, x.slice(x_start, x_size));
}

void SurfaceTarget::update(Vector& x, const SurfaceField& surf) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_state(x.slice(x_start, x_size), surf, type);
}

std::ostream& operator<<(std::ostream& os, const LineTarget&) {
  return os << "Line key value: ";
}

void LineTarget::update(ArrayOfAbsorptionBand& absorption_bands,
                        const Vector& x) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_model(absorption_bands, type, x.slice(x_start, x_size));
}

void LineTarget::update(Vector& x,
                        const ArrayOfAbsorptionBand& absorption_bands) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_state(x.slice(x_start, x_size), absorption_bands, type);
}

const std::vector<AtmTarget>& Targets::atm() const {
  return target<AtmTarget>();
}

const std::vector<SurfaceTarget>& Targets::surf() const {
  return target<SurfaceTarget>();
}

const std::vector<LineTarget>& Targets::line() const {
  return target<LineTarget>();
}

std::vector<AtmTarget>& Targets::atm() { return target<AtmTarget>(); }

std::vector<SurfaceTarget>& Targets::surf() { return target<SurfaceTarget>(); }

std::vector<LineTarget>& Targets::line() { return target<LineTarget>(); }

void Targets::finalize(const AtmField& atmospheric_field,
                       const SurfaceField& surface_field,
                       const ArrayOfAbsorptionBand&) {
  zero_out_x();

  const Size natm  = atm().size();
  const Size nsurf = surf().size();
  const Size nline = line().size();
  ARTS_ASSERT(target_count() == (natm + nsurf + nline));

  Size last_size = 0;

  for (Size i = 0; i < natm; i++) {
    AtmTarget& t = atm()[i];
    ARTS_USER_ERROR_IF(
        std::ranges::any_of(
            atm() | std::views::drop(i + 1), Cmp::eq(t.type), &AtmTarget::type),
        "Multiple targets of the same type: ",
        t.type)
    t.x_start  = last_size;
    t.x_size   = atmospheric_field[t.type].flat_view().size();
    last_size += t.x_size;
  }

  for (Size i = 0; i < nsurf; i++) {
    SurfaceTarget& t = surf()[i];
    ARTS_USER_ERROR_IF(std::ranges::any_of(surf() | std::views::drop(i + 1),
                                           Cmp::eq(t.type),
                                           &SurfaceTarget::type),
                       "Multiple targets of the same type: ",
                       t.type)
    t.x_start  = last_size;
    t.x_size   = surface_field[t.type].flat_view().size();
    last_size += t.x_size;
  }

  for (Size i = 0; i < nline; i++) {
    LineTarget& t = line()[i];
    ARTS_USER_ERROR_IF(std::ranges::any_of(line() | std::views::drop(i + 1),
                                           Cmp::eq(t.type),
                                           &LineTarget::type),
                       "Multiple targets of the same type: ",
                       t.type)
    t.x_start  = last_size;
    t.x_size   = 1;
    last_size += t.x_size;
  }

  finalized = true;

  throwing_check(last_size);
}

std::ostream& operator<<(std::ostream& os, const Targets& targets) {
  os << "Jacobian targets:\n";
  for (const auto& t : targets.atm()) {
    os << "  " << t << '\n';
  }

  for (const auto& t : targets.surf()) {
    os << "  " << t << '\n';
  }

  for (const auto& t : targets.line()) {
    os << "  " << t << '\n';
  }
  return os;
}
}  // namespace Jacobian
