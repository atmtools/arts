#include <workspace.h>

#include <cstdlib>
#include <stdexcept>
#include <tuple>
#include <utility>

int main() try {
  const Numeric*                             returned_fit_data = nullptr;
  const Numeric*                             returned_jac_data = nullptr;
  const JacobianTargets                      targets;
  const Index                                do_jac = 1;
  const measurement_inversion_agendaOperator op{[&returned_fit_data, &returned_jac_data, &targets, &do_jac](
                                                    const JacobianTargets& input_targets, const Index& input_do_jac) {
    if (&input_targets != &targets || &input_do_jac != &do_jac) {
      throw std::runtime_error("Operator input references were copied");
    }
    Vector fit{1, 2, 3};
    Matrix jac(3, 2);
    jac               = 4;
    returned_fit_data = fit.data_handle();
    returned_jac_data = jac.data_handle();
    return std::tuple{std::move(fit), std::move(jac)};
  }};

  // Repeated calls replace populated outputs as well as initially empty ones.
  Vector fit;
  Matrix jac;
  for (Index call = 0; call < 2; ++call) {
    measurement_inversion_agendaExecuteOperator(fit, jac, targets, do_jac, op);
    if (fit.data_handle() != returned_fit_data || jac.data_handle() != returned_jac_data) {
      throw std::runtime_error("Operator outputs were copied instead of transferring their allocations");
    }
    if (fit.size() != 3 || fit[2] != 3 || jac.nrows() != 3 || jac.ncols() != 2 || jac[2, 1] != 4) {
      throw std::runtime_error("Transferred operator outputs have incorrect values");
    }
  }
  return EXIT_SUCCESS;
} catch (const std::exception& error) {
  std::println(stderr, "Agenda operator regression: {}", error.what());
  return EXIT_FAILURE;
}
