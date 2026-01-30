#pragma once

#include <lbl_data.h>
#include <matpack.h>
#include <partfun.h>

#include <Faddeeva.hh>

//! FIXME: These functions should be elsewhere?
namespace Jacobian {
struct Targets;
}  // namespace Jacobian

namespace lbl::voigt::lte::matrix {
void str_scale(ComplexVectorView a,
               const AtmPoint& atm,
               const ConstVectorView& f);

void sumup(ComplexVectorView a,
           const ConstMatrixView& mat,
           const ConstVectorView& f);

void prepare(Matrix& mat,
             const AtmPoint& atm,
             const AbsorptionBands& bands,
             const ZeemanPolarization& pol);

void str_scale(ComplexMatrixView a,
               const AtmPoint& atm,
               const ConstVectorView& fs,
               const std::vector<bool>& df,
               const Size it);

void sumup(ComplexMatrixView a,
           const ConstMatrixView& mat,
           const ConstVectorView& f,
           const std::vector<bool>& df);

void prepare(Matrix& mat,
             const AtmPoint& atm,
             const AbsorptionBands& bands,
             const Jacobian::Targets& jac_targets,
             const ZeemanPolarization& pol);
}  // namespace lbl::voigt::lte::matrix
