#include "quantum_term_symbol.h"

#include "debug.h"

namespace Quantum::Helpers {
namespace detail {
constexpr std::string_view select_molecular_lambda(Index lambda) {
  switch (lambda) {
    case 0:
      return R"(\Sigma)";
    case 1:
      return R"(\Gamma)";
    case 2:
      return R"(\Delta)";
    case 3:
      return R"(\Phi)";
  }
  ARTS_USER_ERROR("Cannot translate lambda ", lambda, " to symbolic value");
}
}  // namespace detail

std::string molecular_term_symbol(const QuantumIdentifier& qid) {
  std::string upp, low;

  if (qid.val.has(QuantumNumberType::ElecStateLabel)) {
    upp += qid.val[QuantumNumberType::ElecStateLabel].str_upp();
    low += qid.val[QuantumNumberType::ElecStateLabel].str_low();
  }

  if (qid.val.has(QuantumNumberType::S)) {
    upp += var_string("$^{", 2 * qid.val[QuantumNumberType::S].upp() + 1, "}$");
    low += var_string("$^{", 2 * qid.val[QuantumNumberType::S].low() + 1, "}$");
  } else {
    upp += "$^{?}$";
    low += "$^{?}$";
  }

  if (qid.val.has(QuantumNumberType::Lambda)) {
    upp += var_string("${",
                      detail::select_molecular_lambda(
                          qid.val[QuantumNumberType::Lambda].upp().toIndex()),
                      "}$");
    low += var_string("${",
                      detail::select_molecular_lambda(
                          qid.val[QuantumNumberType::Lambda].low().toIndex()),
                      "}$");
  } else {
    upp += "$?$";
    low += "$?$";
  }

  const bool vibInv = qid.val.has(QuantumNumberType::vibInv);
  const bool Omega = qid.val.has(QuantumNumberType::Omega);
  if (Omega) {
    upp += var_string("$_{", qid.val[QuantumNumberType::Omega].upp());
    low += var_string("$_{", qid.val[QuantumNumberType::Omega].low());

    if (not vibInv) {
      upp += "}$";
      low += "}$";
    }
  }

  if (vibInv) {
    if (not Omega) {
      upp +=
          var_string("$_{", qid.val[QuantumNumberType::vibInv].str_upp(), "}$");
      low +=
          var_string("$_{", qid.val[QuantumNumberType::vibInv].str_low(), "}$");
    } else {
      upp +=
          var_string(",", qid.val[QuantumNumberType::vibInv].str_upp(), "}$");
      low +=
          var_string(",", qid.val[QuantumNumberType::vibInv].str_low(), "}$");
    }
  }

  if (qid.val.has(QuantumNumberType::parity)) {
    upp +=
        var_string("$^{", qid.val[QuantumNumberType::parity].str_upp(), "}$");
    low +=
        var_string("$^{", qid.val[QuantumNumberType::parity].str_low(), "}$");
  }


  if (qid.val.has(QuantumNumberType::v)) {
    upp +=
        var_string("$\\left(\\nu=", qid.val[QuantumNumberType::v].upp(), "\\right)$");
    low +=
        var_string("$\\left(\\nu=", qid.val[QuantumNumberType::v].low(), "\\right)$");
  }

  return var_string(low, R"($ \leftarrow $)", upp);
}
}  // namespace Quantum::Helpers
