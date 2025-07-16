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

  if (qid.state.contains(QuantumNumberType::ElecStateLabel)) {
    upp += qid.state.at(QuantumNumberType::ElecStateLabel).upper.get<String>();
    low += qid.state.at(QuantumNumberType::ElecStateLabel).lower.get<String>();
  }

  if (qid.state.contains(QuantumNumberType::S)) {
    upp +=
        std::format("$^{{{}}}$", 2 * qid.state.at(QuantumNumberType::S).upper.get<Rational>() + 1);
    low +=
        std::format("$^{{{}}}$", 2 * qid.state.at(QuantumNumberType::S).lower.get<Rational>() + 1);
  } else {
    upp += "$^{?}$";
    low += "$^{?}$";
  }

  if (qid.state.contains(QuantumNumberType::Lambda)) {
    upp += std::format("${{{}}}$",
                       select_molecular_lambda(
                           qid.state.at(QuantumNumberType::Lambda).upper.get<Rational>().toIndex()));
    low += std::format("${{{}}}$",
                       select_molecular_lambda(
                           qid.state.at(QuantumNumberType::Lambda).lower.get<Rational>().toIndex()));
  } else {
    upp += "$?$";
    low += "$?$";
  }

  const bool vibInv = qid.state.contains(QuantumNumberType::vibInv);
  const bool Omega  = qid.state.contains(QuantumNumberType::Omega);
  if (Omega) {
    upp += std::format("$_{{{}", qid.state.at(QuantumNumberType::Omega).upper.get<Rational>());
    low += std::format("$_{{{}", qid.state.at(QuantumNumberType::Omega).lower.get<Rational>());

    if (not vibInv) {
      upp += "}$";
      low += "}$";
    }
  }

  if (vibInv) {
    if (not Omega) {
      upp += std::format("$_{{{}}}$",
                         qid.state.at(QuantumNumberType::vibInv).upper.get<String>());
      low += std::format("$_{{{}}}$",
                         qid.state.at(QuantumNumberType::vibInv).lower.get<String>());
    } else {
      upp +=
          std::format(",{}}}$", qid.state.at(QuantumNumberType::vibInv).upper.get<String>());
      low +=
          std::format(",{}}}$", qid.state.at(QuantumNumberType::vibInv).lower.get<String>());
    }
  }

  if (qid.state.contains(QuantumNumberType::parity)) {
    upp +=
        std::format("$^{{{}}}$", qid.state.at(QuantumNumberType::parity).upper.get<String>());
    low +=
        std::format("$^{{{}}}$", qid.state.at(QuantumNumberType::parity).lower.get<String>());
  }

  if (qid.state.contains(QuantumNumberType::v)) {
    upp += std::format(R"($\left(\nu={}\right)$)",
                       qid.state.at(QuantumNumberType::v).upper.get<Rational>());
    low += std::format(R"($\left(\nu={}\right)$)",
                       qid.state.at(QuantumNumberType::v).lower.get<Rational>());
  }

  return std::format(R"({}$ \leftarrow ${})", low, upp);
}
}  // namespace Quantum::Helpers
