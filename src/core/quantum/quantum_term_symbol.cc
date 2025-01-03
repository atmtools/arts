#include "quantum_term_symbol.h"

#include <debug.h>

namespace {
constexpr std::string_view select_molecular_lambda(Index lambda) {
  switch (lambda) {
    case 0: return R"(\Sigma)";
    case 1: return R"(\Gamma)";
    case 2: return R"(\Delta)";
    case 3: return R"(\Phi)";
  }
  ARTS_USER_ERROR("Cannot translate lambda {} to symbolic value", lambda);
}
}  // namespace

namespace Quantum::Helpers {
std::string molecular_term_symbol(const QuantumIdentifier& qid) {
  std::string upp, low;

  if (qid.val.has(QuantumNumberType::ElecStateLabel)) {
    upp += qid.val[QuantumNumberType::ElecStateLabel].str_upp();
    low += qid.val[QuantumNumberType::ElecStateLabel].str_low();
  }

  if (qid.val.has(QuantumNumberType::S)) {
    upp +=
        std::format("$^{{{}}}$", 2 * qid.val[QuantumNumberType::S].upp() + 1);
    low +=
        std::format("$^{{{}}}$", 2 * qid.val[QuantumNumberType::S].low() + 1);
  } else {
    upp += "$^{?}$";
    low += "$^{?}$";
  }

  if (qid.val.has(QuantumNumberType::Lambda)) {
    upp += std::format("${{{}}}$",
                       select_molecular_lambda(
                           qid.val[QuantumNumberType::Lambda].upp().toIndex()));
    low += std::format("${{{}}}$",
                       select_molecular_lambda(
                           qid.val[QuantumNumberType::Lambda].low().toIndex()));
  } else {
    upp += "$?$";
    low += "$?$";
  }

  const bool vibInv = qid.val.has(QuantumNumberType::vibInv);
  const bool Omega  = qid.val.has(QuantumNumberType::Omega);
  if (Omega) {
    upp += std::format("$_{{{}", qid.val[QuantumNumberType::Omega].upp());
    low += std::format("$_{{{}", qid.val[QuantumNumberType::Omega].low());

    if (not vibInv) {
      upp += "}$";
      low += "}$";
    }
  }

  if (vibInv) {
    if (not Omega) {
      upp += std::format("$_{{{}}}$",
                         qid.val[QuantumNumberType::vibInv].str_upp());
      low += std::format("$_{{{}}}$",
                         qid.val[QuantumNumberType::vibInv].str_low());
    } else {
      upp +=
          std::format(",{}}}$", qid.val[QuantumNumberType::vibInv].str_upp());
      low +=
          std::format(",{}}}$", qid.val[QuantumNumberType::vibInv].str_low());
    }
  }

  if (qid.val.has(QuantumNumberType::parity)) {
    upp +=
        std::format("$^{{{}}}$", qid.val[QuantumNumberType::parity].str_upp());
    low +=
        std::format("$^{{{}}}$", qid.val[QuantumNumberType::parity].str_low());
  }

  if (qid.val.has(QuantumNumberType::v)) {
    upp += std::format(R"($\left(\nu={}\right)$)",
                       qid.val[QuantumNumberType::v].upp());
    low += std::format(R"($\left(\nu={}\right)$)",
                       qid.val[QuantumNumberType::v].low());
  }

  return std::format(R"({}$ \leftarrow ${})", low, upp);
}
}  // namespace Quantum::Helpers
