#pragma once

#include <lbl.h>
#include "fwd_cia.h"
#include "fwd_hxsec.h"
#include "fwd_predef.h"
#include "rtepack.h"

namespace fwd {
struct propmat_operator {
lbl::fwd::line_storage lines{};
cia::full cia{};
predef::full predef{};
hxsec::full hxsec{};

propmat_operator() = default;
propmat_operator(const propmat_operator&) = default;
propmat_operator(propmat_operator&&) = default;
propmat_operator& operator=(const propmat_operator&) = default;
propmat_operator& operator=(propmat_operator&&) = default;

propmat_operator(const lbl::fwd::line_storage& lines,
                 const cia::full& cia,
                 const predef::full& predef,
                 const hxsec::full& hxsec);

std::pair<Propmat, Stokvec> operator()(const Numeric frequency, const Vector2 los={0, 0}) const;
};  // struct propmat_operator
}  // namespace fwd
