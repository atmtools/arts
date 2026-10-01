#include <arts_constants.h>
#include <arts_conversions.h>
#include <atm.h>
#include <jacobian.h>
#include <lbl_lineshape_voigt_ecs.h>
#include <lbl_lineshape_voigt_ecs_hartmann.h>
#include <lbl_lineshape_voigt_ecs_makarov.h>
#include <lbl_lineshape_voigt_lte.h>
#include <physics_funcs.h>
#include <wigner_functions.h>

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {
void require(bool ok, std::string_view message) {
  if (not ok) throw std::runtime_error(std::string(message));
}
void near(Complex actual, Complex expected, Numeric tolerance, std::string_view message) {
  require(std::abs(actual - expected) <= tolerance * std::max(std::abs(expected), Numeric{1e-100}), message);
}
template <class F> void throws(F&& f, std::string_view text) {
  try {
    f();
  } catch (const std::exception& e) {
    require(std::string_view(e.what()).contains(text), e.what());
    return;
  }
  throw std::runtime_error(std::format("Expected an error containing '{}'", text));
}

QuantumIdentifier co2_id() {
  QuantumIdentifier id{"CO2-626"_isot};
  id.state[QuantumNumberType::l2] = {.upper = Rational{0}, .lower = Rational{0}};
  return id;
}

lbl::band_data co2_band() {
  lbl::band_data band;
  band.lineshape = LineByLineLineshape::VP_ECS_HARTMANN;
  lbl::line ln;
  ln.f0                       = 1e11;
  ln.a                        = 1e-5;
  ln.gu                       = 3;
  ln.gl                       = 1;
  ln.ls.T0                    = 296;
  ln.qn[QuantumNumberType::J] = {.upper = Rational{1}, .lower = Rational{0}};
  ln.ls.single_models[SpeciesEnum::Nitrogen].data[LineShapeModelVariable::G0] =
      lbl::temperature::data{LineShapeModelType::T0, {1.0}};
  band.lines.push_back(ln);
  return band;
}

LinemixingSpeciesEcsData collision_data() {
  LinemixingSingleEcsData rate;
  rate.scaling              = {LineShapeModelType::T0, {1.0}};
  rate.beta                 = {LineShapeModelType::T0, {0.01}};
  rate.lambda               = {LineShapeModelType::T0, {0.6}};
  rate.collisional_distance = {LineShapeModelType::T0, {2e-10}};
  return {{SpeciesEnum::Nitrogen, rate}, {SpeciesEnum::Oxygen, rate}};
}

void isolated_line() {
  const auto id    = co2_id();
  const auto band  = co2_band();
  const auto rates = collision_data();
  for (Numeric pressure : {1.0, 1e5, 1e7}) {
    AtmPoint atm;
    atm.pressure                    = pressure;
    atm.temperature                 = 296;
    atm[SpeciesEnum::Nitrogen]      = 1;
    atm[SpeciesEnum::CarbonDioxide] = 4e-4;
    atm[id.isot]                    = 1;
    const Numeric D =
        band.front().f0 * std::sqrt(Constant::doppler_broadening_const_squared * atm.temperature / id.isot.mass);
    Vector                       grid{band.front().f0 - 2 * D, band.front().f0, band.front().f0 + D};
    lbl::voigt::ecs::ComputeData data(grid, atm);
    PropmatVector                pm(grid.size());
    pm = 0;
    PropmatMatrix     dpm(0, grid.size());
    Jacobian::Targets targets;
    lbl::voigt::ecs::calculate(
        pm, dpm, data, grid, Range(0, grid.size()), targets, id, band, rates, atm, ZeemanPolarization::no, false);
    lbl::voigt::lte::single_shape isolated(id.isot, band.front(), atm, ZeemanPolarization::no, 0);
    for (Size i = 0; i < grid.size(); ++i) {
      const Numeric factor = data.scl[i] * Constant::c * Constant::c / (8 * Constant::pi);
      near(pm[i].A(), factor * isolated(grid[i]).real(), 2e-12, "Isolated ECS differs from ordinary Voigt");
    }
    require(data.sum_rule_residual[0] == 1, "Single-line closure residual should expose the missing couplings");
  }
}

void eigen_resolvent() {
  AtmPoint atm;
  atm.pressure    = 1e5;
  atm.temperature = 296;
  Vector                       grid{1e9};
  lbl::voigt::ecs::ComputeData data(grid, atm);
  data.pop  = Vector{2.0, 0.5};
  data.dip  = Vector{1.0, -0.8};
  data.vmrs = Vector{1.0};
  data.Ws.resize(1, 2, 2);
  data.Ws[0][0, 0] = Complex(1e9, 2);
  data.Ws[0][1, 1] = Complex(1e9 + 5, 3);
  data.Ws[0][0, 1] = Complex(0, -0.2);
  data.Ws[0][1, 0] = Complex(0, -0.1);
  const ComplexMatrix original{data.Ws[0]};
  data.core_calc_eqv();
  const ComplexMatrix strengths{data.eqv_strs};
  const ComplexMatrix positions{data.eqv_vals};
  // core_calc_eqv subtracts the carrier from Ws in place, exactly as it does
  // in production, where adapt() rebuilds Ws before every call. Restore the
  // input here so that repeating the calculation repeats the same problem.
  data.Ws[0] = original;
  data.core_calc_eqv();
  Complex total = 0;
  for (Index i = 0; i < 2; ++i) {
    near(data.eqv_strs[0, i], strengths[0, i], 1e-14, "Equivalent strengths changed on repeated calculation");
    near(data.eqv_vals[0, i], positions[0, i], 1e-14, "Equivalent positions changed on repeated calculation");
    total += data.eqv_strs[0, i];
  }
  near(total, 2.0 + 0.5 * 0.8 * 0.8, 2e-14, "Equivalent strengths violate residue-sum invariant");
  require(data.eigenvector_rcond[0] > 0.1, "Unexpectedly ill-conditioned test eigenvectors");
  for (Numeric frequency : {1e9 - 4, 1e9 + 2, 1e9 + 9}) {
    ComplexMatrix A(2, 2);
    ComplexVector rhs(2), solution(2);
    for (Index i = 0; i < 2; ++i) {
      rhs[i] = data.pop[i] * data.dip[i];
      for (Index j = 0; j < 2; ++j) A[i, j] = (i == j ? Complex(frequency) : Complex(0)) - original[i, j];
    }
    solve(solution, A, rhs);
    Complex direct = 0, equivalent = 0;
    for (Index i = 0; i < 2; ++i) {
      direct     += data.dip[i] * solution[i];
      equivalent += data.eqv_strs[0, i] / (frequency - data.eqv_vals[0, i]);
    }
    near(equivalent, direct, 5e-8, "Equivalent-line resolvent differs from direct solve");
  }
  // Each case seeds the whole matrix, since a failed call also leaves Ws centred.
  data.Ws[0][0, 0] = Complex(1e9, 2);
  data.Ws[0][1, 1] = Complex(1e9, 2);
  data.Ws[0][0, 1] = Complex(0, 1);
  data.Ws[0][1, 0] = 0;
  throws([&] { data.core_calc_eqv(); }, "condition");
  data.Ws[0][0, 0] = Complex(1e9, -1);
  data.Ws[0][1, 1] = Complex(1e9, 2);
  data.Ws[0][0, 1] = 0;
  data.Ws[0][1, 0] = 0;
  throws([&] { data.core_calc_eqv(); }, "negative damping");
}

void makarov_and_presorting() {
  QuantumIdentifier id{"O2-66"_isot};
  id.state[QuantumNumberType::S]        = {.upper = Rational{1}, .lower = Rational{1}};
  auto band                             = co2_band();
  band.lineshape                        = LineByLineLineshape::VP_ECS_MAKAROV;
  band.front().qn[QuantumNumberType::N] = {.upper = Rational{1}, .lower = Rational{1}};
  AtmPoint atm;
  atm.pressure               = 1e5;
  atm.temperature            = 296;
  atm[SpeciesEnum::Nitrogen] = 1;
  Vector                       grid{band.front().f0};
  lbl::voigt::ecs::ComputeData data(grid, atm);
  data.adapt_single(id, band, collision_data(), atm);
  data.core_calc(grid);
  require(data.shape[0].real() > 0, "Supported Makarov single line has no absorption");
  atm.specs.clear();
  id.state[QuantumNumberType::S].upper = Rational{0};
  throws([&] { lbl::voigt::ecs::validate_band(id, band); }, "S=1");
  id.state[QuantumNumberType::S].upper        = Rational{1};
  band.front().qn[QuantumNumberType::N].lower = Rational{3};
  throws([&] { lbl::voigt::ecs::validate_band(id, band); }, "N");

  band = co2_band();
  for (Index i = 1; i < 3; ++i) {
    auto ln                      = band.front();
    ln.f0                       += Numeric(i) * 1e7;
    ln.a                        *= Numeric(i + 1);
    ln.qn[QuantumNumberType::J]  = {.upper = Rational{2 * i + 1}, .lower = Rational{2 * i}};
    band.lines.push_back(ln);
  }
  const auto co2 = co2_id();
  Vector     expected_pop(3), expected_dip(3), expected_dipr(3);
  for (Index i = 0; i < 3; ++i) {
    lbl::band_data single;
    single.lineshape = band.lineshape;
    single.lines.push_back(band.lines[i]);
    data.adapt_single(co2, single, {}, atm);
    expected_pop[i]  = data.pop[0];
    expected_dip[i]  = data.dip[0];
    expected_dipr[i] = data.dipr[0];
  }
  data.sort = ArrayOfIndex{1, 2, 0};
  data.adapt_single(co2, band, {}, atm, true);
  for (Index i = 0; i < 3; ++i) {
    const Index j = data.sort[i];
    near(data.pop[i], expected_pop[j], 1e-14, "Presorted population mismatch");
    near(data.dip[i], expected_dip[j], 1e-14, "Presorted dipole mismatch");
    near(data.dipr[i], expected_dipr[j], 1e-14, "Presorted reduced dipole mismatch");
    near(data.Ws[0][i, i].real(), band.lines[j].f0, 1e-14, "Presorted line position mismatch");
  }
}

void validation_and_mixtures() {
  auto       band  = co2_band();
  const auto id    = co2_id();
  const auto rates = collision_data();
  AtmPoint   atm;
  atm.pressure               = 1e5;
  atm.temperature            = 296;
  atm[SpeciesEnum::Nitrogen] = 0.8;
  atm[SpeciesEnum::Oxygen]   = 0.2;
  Vector                       grid{1e11};
  lbl::voigt::ecs::ComputeData data(grid, atm);
  auto                         oxygen                = band.front().ls.single_models.at(SpeciesEnum::Nitrogen);
  oxygen.data[LineShapeModelVariable::G0]            = {LineShapeModelType::T0, {3.0}};
  oxygen.data[LineShapeModelVariable::D0]            = {LineShapeModelType::T0, {0.1}};
  band.front().ls.single_models[SpeciesEnum::Oxygen] = oxygen;
  data.adapt_single(id, band, rates, atm);
  const Complex mixed = data.Ws[0][0, 0];
  data.adapt_multi(id, band, rates, atm);
  Complex weighted = 0;
  for (Size i = 0; i < data.vmrs.size(); ++i) weighted += data.vmrs[i] * data.Ws[i][0, 0];
  near(mixed, weighted, 1e-14, "Single ECS mixture differs from sum of partner matrices");
  near(mixed.imag(), 1.4e5, 1e-14, "Wrong mixture broadening");

  // Unordered maps with the same keys need not have the same iteration order.
  auto next                      = band.front();
  next.f0                       += 1e7;
  next.qn[QuantumNumberType::J]  = {.upper = Rational{3}, .lower = Rational{2}};
  next.ls.single_models.clear();
  next.ls.single_models[SpeciesEnum::Oxygen]   = oxygen;
  next.ls.single_models[SpeciesEnum::Nitrogen] = band.front().ls.single_models.at(SpeciesEnum::Nitrogen);
  band.lines.push_back(next);
  data.adapt_single(id, band, rates, atm);
  band.lines.back().ls.single_models.erase(SpeciesEnum::Oxygen);
  throws([&] { lbl::voigt::ecs::validate_band(id, band); }, "same broadening species");

  band = co2_band();
  band.lines.push_back(band.front());
  throws([&] { lbl::voigt::ecs::validate_band(id, band); }, "same rotational pair");
  band = co2_band();
  atm.specs.clear();
  auto unsupported_id = id;
  unsupported_id.isot = "CO2-628"_isot;
  throws([&] { data.adapt_single(unsupported_id, band, {}, atm); }, "CO2-626");
  data.adapt_single(id, band, {}, atm);
  require(data.Ws[0][0, 0].imag() == 0, "Absent perturbers should produce zero width");
  data.core_calc(grid);
  require(std::isfinite(data.shape[0].real()), "Absent perturbers produced non-finite spectrum");
  atm[SpeciesEnum::Nitrogen] = -0.1;
  throws([&] { data.adapt_single(id, band, rates, atm); }, "VMR");

  ComplexTensor3 strength(1, 1, 1), position(2, 1, 1);
  Vector         temperatures{296};
  throws([&] { lbl::voigt::ecs::equivalent_values(strength, position, data, id, band, rates, atm, temperatures); },
         "same shape");
}

void explicit_partner_mixture() {
  const auto id   = co2_id();
  auto       band = co2_band();
  for (Index i = 1; i < 3; ++i) {
    auto ln                      = band.front();
    ln.f0                       += Numeric(i) * 1e7;
    ln.a                        *= Numeric(i + 1);
    ln.qn[QuantumNumberType::J]  = {.upper = Rational{2 * i + 1}, .lower = Rational{2 * i}};
    band.lines.push_back(ln);
  }
  for (Size i = 0; i < band.size(); ++i) {
    auto& models                            = band.lines[i].ls.single_models;
    auto  oxygen                            = models.at(SpeciesEnum::Nitrogen);
    oxygen.data[LineShapeModelVariable::G0] = {LineShapeModelType::T0, {2.0 + Numeric(i)}};
    oxygen.data[LineShapeModelVariable::D0] = {LineShapeModelType::T0, {0.1 * Numeric(i + 1)}};
    models[SpeciesEnum::Oxygen]             = oxygen;
  }
  auto rates                                         = collision_data();
  rates.at(SpeciesEnum::Oxygen).beta                 = {LineShapeModelType::T0, {0.05}};
  rates.at(SpeciesEnum::Oxygen).lambda               = {LineShapeModelType::T0, {0.8}};
  rates.at(SpeciesEnum::Oxygen).collisional_distance = {LineShapeModelType::T0, {3e-10}};
  AtmPoint atm;
  atm.pressure               = 1e5;
  atm[SpeciesEnum::Nitrogen] = 0.79;
  atm[SpeciesEnum::Oxygen]   = 0.21;
  Vector grid{band.front().f0};
  for (Numeric temperature : {220.0, 296.0, 330.0}) {
    atm.temperature = temperature;
    lbl::voigt::ecs::ComputeData data(grid, atm);
    data.adapt_single(id, band, rates, atm);
    const ComplexMatrix mixed{data.Ws[0]};
    data.adapt_multi(id, band, rates, atm);
    bool has_coupling = false;
    for (Size r = 0; r < band.size(); ++r) {
      for (Size c = 0; c < band.size(); ++c) {
        Complex expected = 0;
        for (Size p = 0; p < data.vmrs.size(); ++p) expected += data.vmrs[p] * data.Ws[p][r, c];
        near(mixed[r, c], expected, 2e-14, "ECS mixture must sum the full collision-partner matrices");
        has_coupling |= r != c and std::abs(mixed[r, c]) > 0;
      }
    }
    require(has_coupling, "Mixture fixture must exercise off-diagonal ECS couplings");
  }
}

void oxygen_energy_levels() {
  using namespace lbl::voigt::ecs;
  require(makarov::level_energy(Rational{1}, Rational{0}) == 0, "O2 ground-state energy must be zero");

  // Independent spectroscopic reference: measured 1- at 118.750340 GHz;
  // the other centers are from the Tretyakov et al. (2005) table used in TRE05.cc
  // (doi:10.1016/j.jms.2004.11.011). The retained approximate level model is
  // accurate only to about 25 MHz for these other branches, not catalog accuracy.
  struct Reference {
    Index   N;
    Numeric minus, plus;
  };
  constexpr std::array reference{Reference{1, 118.750340e9, 56.264774e9},
                                 Reference{3, 62.486253e9, 58.446588e9},
                                 Reference{5, 60.306056e9, 59.590983e9}};
  for (const auto& [n, minus, plus] : reference) {
    const Rational N{n};
    const Numeric  middle      = makarov::level_energy(N, N);
    const Numeric  lower_minus = makarov::level_energy(N, N - 1);
    const Numeric  lower_plus  = makarov::level_energy(N, N + 1);
    const Numeric  tolerance   = n == 1 ? 1e4 : 3e7;
    require(std::abs(Conversion::joule2hz(middle - lower_minus) - minus) < tolerance,
            "O2 N- level splitting disagrees with the spectroscopic reference");
    require(std::abs(Conversion::joule2hz(middle - lower_plus) - plus) < 3e7,
            "O2 N+ level splitting exceeds the approximate model's 30 MHz tolerance");
    near(makarov::rotational_energy(N),
         middle,
         2e-14,
         "O2 reference rotor and resolved levels must use the same energy zero");
    require(lower_minus != lower_plus, "O2 resolved energies must distinguish the two fine-structure branches");
  }

  QuantumIdentifier id{"O2-66"_isot};
  id.state[QuantumNumberType::S] = {.upper = Rational{1}, .lower = Rational{1}};
  std::array<rotational_line, 2> lines{
      {{Rational{3}, Rational{2}, Rational{3}, Rational{3}}, {Rational{3}, Rational{4}, Rational{3}, Rational{3}}}};
  energy_data energies;
  makarov::prepare_energies(energies, id, lines);
  near(energies.e0[0],
       makarov::level_energy(Rational{3}, Rational{2}),
       2e-14,
       "O2 N- preparation lost the resolved lower J");
  near(energies.e0[1],
       makarov::level_energy(Rational{3}, Rational{4}),
       2e-14,
       "O2 N+ preparation lost the resolved lower J");
  require(energies.e0[0] != energies.e0[1], "O2 preparation collapsed distinct fine-structure branches");
  std::swap(lines[0], lines[1]);
  energy_data reversed;
  makarov::prepare_energies(reversed, id, lines);
  near(reversed.e0[0], energies.e0[1], 2e-14, "O2 resolved energies must follow line ordering");
  near(reversed.e0[1], energies.e0[0], 2e-14, "O2 resolved energies must follow line ordering");
  require(reversed.rotational.size() == energies.rotational.size(), "Sorting changed the O2 basis size");
  for (Size i = 0; i < energies.rotational.size(); ++i) {
    near(reversed.rotational[i], energies.rotational[i], 2e-14, "Sorting changed the O2 reference ladder");
    near(reversed.rotational_minus_two[i],
         energies.rotational_minus_two[i],
         2e-14,
         "Sorting changed the O2 shifted reference ladder");
  }
}

void sum_rule_energy_preparation() {
  AtmPoint atm;
  atm.pressure               = 1e5;
  atm.temperature            = 296;
  atm[SpeciesEnum::Nitrogen] = 1;
  const auto   rates         = collision_data();
  const Vector catalogue_e0{Constant::k * 100, Constant::k * 150, Constant::k * 125};
  const Vector rotational_e0{0, Conversion::kaycm2joule(0.39021 * 6), Conversion::kaycm2joule(0.39021 * 20)};

  for (Index model = 0; model < 3; ++model) {
    const bool makarov = model == 2;
    auto       id      = makarov ? QuantumIdentifier{"O2-66"_isot} : co2_id();
    auto       band    = co2_band();
    if (makarov) {
      id.state[QuantumNumberType::S] = {.upper = Rational{1}, .lower = Rational{1}};
      band.lineshape                 = LineByLineLineshape::VP_ECS_MAKAROV;
    } else if (model == 1) {
      // The angular kernel swaps its J roles; closure must still use original lower J.
      id.state[QuantumNumberType::l2].upper = Rational{1};
    }
    const auto prototype = band.front();
    band.lines.clear();
    for (Index i = 0; i < 3; ++i) {
      auto ln                      = prototype;
      ln.f0                       += Numeric(i) * 1e7;
      ln.a                        *= Numeric(i + 1);
      ln.e0                        = catalogue_e0[i];
      ln.qn[QuantumNumberType::J]  = {.upper = Rational{2 * i + 1}, .lower = Rational{2 * i}};
      if (makarov) ln.qn[QuantumNumberType::N] = {.upper = Rational{2 * i + 1}, .lower = Rational{2 * i + 1}};
      band.lines.push_back(ln);
    }

    lbl::voigt::ecs::ComputeData data({}, atm);
    const auto                   check = [&] {
      bool has_coupling = false;
      for (Index i = 0; i < 3; ++i) {
        const auto&   line = band.lines[data.sort[i]];
        const auto&   qn   = data.rotational_lines[i];
        const Numeric expected =
            makarov ? lbl::voigt::ecs::makarov::level_energy(qn.Nl, qn.Jl) : rotational_e0[data.sort[i]];
        near(data.energies.e0[i], expected, 2e-14, "Wrong prepared sum-rule energy or ordering");
        near(data.Ws[0][i, i].real(), line.f0, 2e-14, "Energy preparation changed the catalogue line frequency");
        const auto& J = line.qn.at(QuantumNumberType::J);
        require(qn.Ju == static_cast<Rational>(J.upper) and qn.Jl == static_cast<Rational>(J.lower),
                "Prepared J does not follow matrix ordering");
        if (makarov) {
          const auto& N = line.qn.at(QuantumNumberType::N);
          require(qn.Nu == static_cast<Rational>(N.upper) and qn.Nl == static_cast<Rational>(N.lower),
                  "Prepared N does not follow matrix ordering");
        }
        require(band.lines[i].e0 == catalogue_e0[i], "Preparing sum-rule energies changed catalogue energies");
        for (Index j = i + 1; j < 3; ++j) {
          const Numeric reverse  = data.Ws[0][i, j].imag();
          has_coupling          |= reverse != 0;
          const Numeric balance =
              std::exp((data.energies.e0[i] - data.energies.e0[j]) / (Constant::k * atm.temperature));
          near(data.Ws[0][j, i].imag(), reverse * balance, 2e-14, "Sum-rule detailed balance uses wrong energies");
        }
      }
      require(has_coupling, "Sum-rule energy fixture needs coupled lines");
      const auto& energies = data.energies;
      require(energies.e0.size() == band.size() and
                  energies.rotational.size() == energies.rotational_minus_two.size() and energies.rotational.size() > 5,
              "Prepared energy dimensions do not cover the line and angular bases");
      for (Size i = 0; i < energies.rotational.size(); ++i) {
        const Rational L{static_cast<Index>(i)};
        const auto     rotor =
            makarov ? lbl::voigt::ecs::makarov::rotational_energy : lbl::voigt::ecs::hartmann::rotational_energy;
        near(energies.rotational[i], rotor(L), 2e-14, "Reference-rotor energy must be indexed by angular momentum");
        near(energies.rotational_minus_two[i],
             rotor(L - 2),
             2e-14,
             "Shifted reference-rotor energy must use angular momentum minus two");
      }
    };
    data.adapt_single(id, band, rates, atm);
    require(data.sort[0] != 0, "Sum-rule energy fixture needs nontrivial strength sorting");
    check();

    const Vector        original_pop{data.pop};
    const ComplexMatrix original_matrix{data.Ws[0]};
    const ArrayOfIndex  original_sort{data.sort};
    auto                shifted = band;
    const Numeric       offset  = Constant::k * 83;
    for (auto& ln : shifted.lines) ln.e0 += offset;
    data.adapt_single(id, shifted, rates, atm);
    for (Index i = 0; i < 3; ++i) {
      require(data.sort[i] == original_sort[i], "Common energy offset changed line ordering");
      near(data.pop[i],
           original_pop[i] * std::exp(-offset / (Constant::k * atm.temperature)),
           2e-14,
           "Optical populations must retain absolute catalogue energies");
      for (Index j = 0; j < 3; ++j)
        near(data.Ws[0][i, j], original_matrix[i, j], 5e-14, "Common energy offset changed collision matrix");
    }

    data.sort = ArrayOfIndex{1, 2, 0};
    data.adapt_single(id, band, rates, atm, true);
    check();

    const auto kernel = makarov ? lbl::voigt::ecs::makarov::relaxation_matrix_offdiagonal
                                : lbl::voigt::ecs::hartmann::relaxation_matrix_offdiagonal;
    const auto matrix = [&](const lbl::voigt::ecs::energy_data& energies) {
      Matrix W(3, 3);
      W = 0;
      for (Index i = 0; i < 3; ++i) W[i, i] = atm.pressure;
      kernel(W,
             id,
             data.rotational_lines,
             band.front().ls.T0,
             SpeciesEnum::Nitrogen,
             rates.at(SpeciesEnum::Nitrogen),
             data.dipr,
             energies,
             atm,
             {},
             {},
             {},
             {});
      return W;
    };
    const Matrix original          = matrix(data.energies);
    auto         changed_energies  = data.energies;
    changed_energies.e0[0]        += offset;
    const Matrix changed           = matrix(changed_energies);
    require(original[0, 1] != 0 and original[0, 2] != 0 and changed[0, 2] != 0,
            "Prepared-energy regression needs two couplings in its first row");
    // Row normalization preserves this raw ratio, exposing the energy input
    // even for Hartmann's swapped angular states.
    near(changed[0, 1] / changed[0, 2],
         original[0, 1] / original[0, 2] * std::exp(offset / (Constant::k * atm.temperature)),
         5e-14,
         "Raw relaxation matrix must use the prepared energies");
    for (Index i = 0; i < 3; ++i) {
      for (Index j = i + 1; j < 3; ++j)
        near(changed[j, i],
             changed[i, j] *
                 std::exp((changed_energies.e0[i] - changed_energies.e0[j]) / (Constant::k * atm.temperature)),
             2e-14,
             "Sum-rule correction must use the same prepared energies as the raw kernel");
    }
  }
}

void undefined_basis_rate() {
  // Q(0) is undefined. Instead of testing for it on every evaluation, check the
  // invariant that replaces it: prepare_basis never asks for it.
  AtmPoint atm;
  atm.temperature            = 296;
  atm.pressure               = 1e5;
  atm[SpeciesEnum::Nitrogen] = 1;
  lbl::voigt::ecs::energy_data energies;
  lbl::voigt::ecs::prepare_rotational_ladder(energies, 4, lbl::voigt::ecs::hartmann::rotational_energy);
  const auto basis = lbl::voigt::ecs::prepare_basis(
      4, energies, collision_data().at(SpeciesEnum::Nitrogen), 296, "CO2-626"_isot, SpeciesEnum::Nitrogen, atm);
  require(basis.Q[0] == 0, "prepare_basis must leave the undefined Q(0) at zero");
  for (Index L = 1; L < 4; ++L) require(basis.Q[L] > 0, "prepare_basis skipped a defined basis rate");
}
}  // namespace

int main() try {
  undefined_basis_rate();
  lbl::voigt::ecs::validate_band(co2_id(), co2_band());
  WignerInformation wigner(100, 0, true, true);
  isolated_line();
  eigen_resolvent();
  makarov_and_presorting();
  validation_and_mixtures();
  explicit_partner_mixture();
  oxygen_energy_levels();
  sum_rule_energy_preparation();
  wigner.unload();
  std::cout << "ECS numerical regression tests passed\n";
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
