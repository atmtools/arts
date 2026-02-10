#include <arts_omp.h>
#include <jacobian.h>
#include <lbl.h>
#include <omp.h>
#include <rng.h>
#include <time_report.h>

#include <iostream>

namespace {
const QuantumIdentifier bnd_qid{"O2-66"_isot};

JacobianTargets create_jac_targets() {
  JacobianTargets jac_targets;

  jac_targets.emplace_back(AtmKey::t);
  jac_targets.emplace_back(AtmKey::p);
  jac_targets.emplace_back(AtmKey::mag_u);
  jac_targets.emplace_back(AtmKey::mag_v);
  jac_targets.emplace_back(AtmKey::mag_w);
  jac_targets.emplace_back(AtmKey::wind_u);
  jac_targets.emplace_back(AtmKey::wind_v);
  jac_targets.emplace_back(AtmKey::wind_w);
  jac_targets.emplace_back(bnd_qid.isot);
  jac_targets.emplace_back(bnd_qid.isot.spec);
  return jac_targets;
}

std::vector<lbl::line> create_lines(Index M) {
  constexpr Index n = 5;
  const Matrix x    = random_numbers<2>({M, n}, 0.0, 1.0);

  std::vector<lbl::line> lines(M);

  for (Index i = 0; i < M; i++) {
    auto& line = lines[i];

    line.a      = 4.479289583303983e-09 * (1 + nonstd::abs(x[i, 0]));
    line.f0     = 118750348044.712 + nonstd::abs(x[i, 1]) * 1e9;
    line.e0     = 1e-23 * (1 + nonstd::abs(x[i, 2]));
    line.gu     = 3.0;
    line.gl     = 1.0;
    line.z.gu() = 1.0011;
    line.z.gl() = 0.0;
    line.ls.T0  = 296.0;
    auto& ls    = line.ls.single_models[bnd_qid.isot.spec];
    ls.data[LineShapeModelVariable::G0] = lbl::temperature::data(
        LineShapeModelType::T1, Vector{x[i, Range(3, 2)]});
    line.qn[QuantumNumberType::J] = {.upper = Rational(1),
                                     .lower = Rational(0)};
  }
  return lines;
}

Numeric lbl_temperature_t0(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_t0");

  using namespace lbl::temperature;
  // constexpr Numeric T0 = 296;
  // constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += model::T0(x[i, 0]);

  return result;
}

Numeric lbl_temperature_t1(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_t1");

  using namespace lbl::temperature;
  constexpr Numeric T0 = 296;
  constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += model::T1(x[i, 0], x[i, 1], T0, T);

  return result;
}

Numeric lbl_temperature_t2(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_t2");

  using namespace lbl::temperature;
  constexpr Numeric T0 = 296;
  constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i)
    result += model::T2(x[i, 0], x[i, 1], x[i, 2], T0, T);

  return result;
}

Numeric lbl_temperature_t3(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_t3");

  using namespace lbl::temperature;
  constexpr Numeric T0 = 296;
  constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += model::T3(x[i, 0], x[i, 1], T0, T);
  return result;
}

Numeric lbl_temperature_t4(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_t4");

  using namespace lbl::temperature;
  constexpr Numeric T0 = 296;
  constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i)
    result += model::T4(x[i, 0], x[i, 1], x[i, 2], T0, T);
  return result;
}

Numeric lbl_temperature_t5(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_t5");

  using namespace lbl::temperature;
  constexpr Numeric T0 = 296;
  constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += model::T5(x[i, 0], x[i, 1], T0, T);
  return result;
}

Numeric lbl_temperature_dpl(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_dpl");

  using namespace lbl::temperature;
  constexpr Numeric T0 = 296;
  constexpr Numeric T  = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i)
    result += model::DPL(x[i, 0], x[i, 1], x[i, 2], x[i, 3], T0, T);
  return result;
}

Numeric lbl_temperature_aer(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_aer");

  using namespace lbl::temperature;
  // constexpr Numeric T0 = 296;
  constexpr Numeric T = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i)
    result += model::AER(x[i, 0], x[i, 1], x[i, 2], x[i, 3], T);
  return result;
}

Numeric lbl_temperature_poly(const Matrix& x) {
  ARTS_NAMED_TIME_REPORT("lbl_temperature_poly");

  using namespace lbl::temperature;
  // constexpr Numeric T0 = 296;
  constexpr Numeric T = 250;

  const Index n = x.nrows();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += model::POLY(x[i], T);
  return result;
}

Numeric lbl_data_line_s(const std::vector<lbl::line>& lines) {
  ARTS_NAMED_TIME_REPORT("lbl_data_line_s");

  const Index n = lines.size();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += lines[i].s(250.0, 182.0);
  return result;
}

Numeric lbl_data_line_ds_da(const std::vector<lbl::line>& lines) {
  ARTS_NAMED_TIME_REPORT("lbl_data_line_ds_da");

  const Index n = lines.size();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += lines[i].ds_da(250.0, 182.0);
  return result;
}

Numeric lbl_data_line_ds_df0(const std::vector<lbl::line>& lines) {
  ARTS_NAMED_TIME_REPORT("lbl_data_line_ds_df0");

  const Index n = lines.size();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += lines[i].ds_df0(250.0, 182.0);
  return result;
}

Numeric lbl_data_line_ds_de0(const std::vector<lbl::line>& lines) {
  ARTS_NAMED_TIME_REPORT("lbl_data_line_ds_de0");

  const Index n = lines.size();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += lines[i].ds_de0(250.0, 182.0);
  return result;
}

Numeric lbl_data_line_ds_dT(const std::vector<lbl::line>& lines) {
  ARTS_NAMED_TIME_REPORT("lbl_data_line_ds_dT");

  const Index n = lines.size();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i) result += lines[i].ds_dT(250.0, 182.0, 72.0);
  return result;
}

Numeric lbl_data_line_hitran_s(const std::vector<lbl::line>& lines) {
  ARTS_NAMED_TIME_REPORT("lbl_data_line_hitran_s");

  const Index n = lines.size();

  Numeric result = 0.0;
  for (Index i = 0; i < n; ++i)
    result += lines[i].hitran_s(bnd_qid.isot, 296.0);
  return result;
}

Numeric lbl_voigt_lte_calculate(PropmatVector& pm,
                                PropmatMatrix& dpm,
                                const AbsorptionBand& bnd,
                                const Vector& fs,
                                const AtmPoint& atm) {
  ARTS_NAMED_TIME_REPORT(std::format("lbl_voigt_lte_calculate; threads: {}",
                                     arts_omp_get_max_threads()));

  using namespace lbl::voigt::lte;

  const ZeemanPolarization pol = ZeemanPolarization::no;
  constexpr Vector2 los{};
  const Jacobian::Targets jac_targets{};

  ComputeData com_data(fs, atm, los, pol);

  const Size n = arts_omp_get_max_threads();
  if (n == 1) {
    calculate(pm,
              dpm,
              com_data,
              fs,
              Range{0, fs.size()},
              jac_targets,
              bnd_qid,
              bnd,
              atm,
              pol,
              false);
  } else {
    const auto f_ranges = matpack::omp_offset_count(fs.size(), n);
#pragma omp parallel for firstprivate(com_data)
    for (Size i = 0; i < n; i++) {
      calculate(pm,
                dpm,
                com_data,
                fs,
                f_ranges[i],
                jac_targets,
                bnd_qid,
                bnd,
                atm,
                pol,
                false);
    }
  }

  return pm[0].A();
}

Numeric lbl_voigt_lte_mirror_calculate(PropmatVector& pm,
                                       PropmatMatrix& dpm,
                                       const AbsorptionBand& bnd,
                                       const Vector& fs,
                                       const AtmPoint& atm) {
  ARTS_NAMED_TIME_REPORT(
      std::format("lbl_voigt_lte_mirror_calculate; threads: {}",
                  arts_omp_get_max_threads()));

  using namespace lbl::voigt::lte_mirror;

  const ZeemanPolarization pol = ZeemanPolarization::no;
  constexpr Vector2 los{};
  const Jacobian::Targets jac_targets{};

  ComputeData com_data(fs, atm, los, pol);

  const Size n = arts_omp_get_max_threads();
  if (n == 1) {
    calculate(pm,
              dpm,
              com_data,
              fs,
              Range{0, fs.size()},
              jac_targets,
              bnd_qid,
              bnd,
              atm,
              pol,
              false);
  } else {
    const auto f_ranges = matpack::omp_offset_count(fs.size(), n);
#pragma omp parallel for firstprivate(com_data)
    for (Size i = 0; i < n; i++) {
      calculate(pm,
                dpm,
                com_data,
                fs,
                f_ranges[i],
                jac_targets,
                bnd_qid,
                bnd,
                atm,
                pol,
                false);
    }
  }

  return pm[0].A();
}

Numeric lbl_voigt_lte_matrix_prepare_manylines(
    Matrix& mat,
    const std::vector<lbl::flat_band_data>& flat,
    const AtmPoint& atm) {
  ARTS_NAMED_TIME_REPORT(
      std::format("lbl_voigt_lte_matrix_prepare_manylines; threads: {}",
                  arts_omp_get_max_threads()));

  lbl::voigt::lte::matrix::prepare(mat, atm, flat, ZeemanPolarization::no);

  return mat[0, 0];
}

Numeric lbl_voigt_lte_matrix_prepare_manylines_jac(
    Matrix& mat,
    const std::vector<lbl::flat_band_data>& flat,
    const AtmPoint& atm,
    const JacobianTargets& jac_targets) {
  ARTS_NAMED_TIME_REPORT(
      std::format("lbl_voigt_lte_matrix_prepare_manylines_jac; threads: {}",
                  arts_omp_get_max_threads()));

  lbl::voigt::lte::matrix::prepare(
      mat, atm, flat, jac_targets, ZeemanPolarization::no);

  return mat[0, 0];
}

Numeric lbl_voigt_lte_matrix_prepare_manybands(
    Matrix& mat,
    const std::vector<lbl::flat_band_data>& flat,
    const AtmPoint& atm) {
  ARTS_NAMED_TIME_REPORT(
      std::format("lbl_voigt_lte_matrix_prepare_manybands; threads: {}",
                  arts_omp_get_max_threads()));

  lbl::voigt::lte::matrix::prepare(mat, atm, flat, ZeemanPolarization::no);

  return mat[0, 0];
}

Numeric lbl_voigt_lte_matrix_prepare_manybands_jac(
    Matrix& mat,
    const std::vector<lbl::flat_band_data>& flat,
    const AtmPoint& atm,
    const JacobianTargets& jac_targets) {
  ARTS_NAMED_TIME_REPORT(
      std::format("lbl_voigt_lte_matrix_prepare_manybands_jac; threads: {}",
                  arts_omp_get_max_threads()));

  lbl::voigt::lte::matrix::prepare(
      mat, atm, flat, jac_targets, ZeemanPolarization::no);

  return mat[0, 0];
}

Numeric lbl_voigt_lte_matrix_prepare_sort(Matrix& mat) {
  ARTS_NAMED_TIME_REPORT(std::format(
      "lbl_voigt_lte_matrix_prepare_sort; mat-shape: {:B,}", mat.shape()));

  lbl::voigt::lte::matrix::sort(mat);

  return mat[0, 0];
}

Numeric lbl_voigt_lte_matrix_sumup(ComplexVectorView a,
                                   const ConstMatrixView mat,
                                   const ConstVectorView f) {
  ARTS_NAMED_TIME_REPORT(std::format("lbl_voigt_lte_matrix_sumup; threads: {}",
                                     arts_omp_get_max_threads()));

  lbl::voigt::lte::matrix::sumup(a, mat, f);

  return a[0].real();
}

Numeric lbl_voigt_lte_matrix_sumup_inf_cutoff(ComplexVectorView a,
                                              const ConstMatrixView mat,
                                              const ConstVectorView f) {
  ARTS_NAMED_TIME_REPORT(
      std::format("lbl_voigt_lte_matrix_sumup_inf_cutoff; threads: {}",
                  arts_omp_get_max_threads()));

  lbl::voigt::lte::matrix::sumup(
      a, mat, f, std::numeric_limits<Numeric>::infinity());

  return a[0].real();
}

Numeric lbl_voigt_lte_matrix_sumup_jac(ComplexMatrixView a,
                                       const ConstMatrixView mat,
                                       const ConstVectorView f,
                                       const std::vector<bool>& df) {
  ARTS_NAMED_TIME_REPORT(
      std::format("lbl_voigt_lte_matrix_sumup_jac; threads: {}",
                  arts_omp_get_max_threads()));

  lbl::voigt::lte::matrix::sumup(a, mat, f, df);

  return a[0, 0].real();
}

Numeric lbl_voigt_lte_matrix_sumup_inf_cutoff_jac(ComplexMatrixView a,
                                                  const ConstMatrixView mat,
                                                  const ConstVectorView f,
                                                  const std::vector<bool>& df) {
  ARTS_NAMED_TIME_REPORT(
      std::format("lbl_voigt_lte_matrix_sumup_inf_cutoff_jac threads: {}",
                  arts_omp_get_max_threads()));

  lbl::voigt::lte::matrix::sumup(
      a, mat, f, std::numeric_limits<Numeric>::infinity(), df);

  return a[0, 0].real();
}
}  // namespace

int main() {
  Numeric buf{};
  {
    constexpr Index M = 10'000'000;
    constexpr Index n = 4;
    const Matrix x    = random_numbers<2>({M, n}, 0.0, 1.0);

    buf += lbl_temperature_t0(x);
    buf += lbl_temperature_t1(x);
    buf += lbl_temperature_t2(x);
    buf += lbl_temperature_t3(x);
    buf += lbl_temperature_t4(x);
    buf += lbl_temperature_t5(x);
    buf += lbl_temperature_dpl(x);
    buf += lbl_temperature_aer(x);
    buf += lbl_temperature_poly(x);
  }

  {
    constexpr Index M                  = 10'000'000;
    const std::vector<lbl::line> lines = create_lines(M);

    buf += lbl_data_line_s(lines);
    buf += lbl_data_line_ds_da(lines);
    buf += lbl_data_line_ds_df0(lines);
    buf += lbl_data_line_ds_de0(lines);
    buf += lbl_data_line_ds_dT(lines);
    buf += lbl_data_line_hitran_s(lines);
  }

  {
    AtmPoint atm;
    constexpr Index M        = 10'000;
    constexpr Index N        = 2'000;
    const AbsorptionBand bnd = {.lines = create_lines(M)};
    atm.temperature          = 250.0;
    atm.pressure             = 182.0;
    atm[bnd_qid.isot.spec]   = 0.21;
    const Vector f           = random_numbers<1>({N}, 0.0, 1.0);
    PropmatVector pm(N);
    PropmatMatrix dpm(0, N);

    buf += lbl_voigt_lte_calculate(pm, dpm, bnd, f, atm);
    buf += lbl_voigt_lte_mirror_calculate(pm, dpm, bnd, f, atm);

    const int cores = arts_omp_get_max_threads();
    arts_omp_set_num_threads(1);

    buf += lbl_voigt_lte_calculate(pm, dpm, bnd, f, atm);
    buf += lbl_voigt_lte_mirror_calculate(pm, dpm, bnd, f, atm);

    arts_omp_set_num_threads(cores);
  }

  {
    constexpr Index M = 1'000'000;
    AbsorptionBands bands{};
    bands[bnd_qid] = AbsorptionBand{.lines = create_lines(M)};
    AtmPoint atm;
    atm.temperature        = 250.0;
    atm.pressure           = 182.0;
    atm[bnd_qid.isot.spec] = 0.21;
    Matrix mat(M, 5);
    const std::vector<lbl::flat_band_data> flat = lbl::flatter_view(
        bands, ZeemanPolarization::no, [](auto&, auto&) { return true; });

    constexpr Index n = 5;
    const Matrix x    = random_numbers<2>({M, n}, 0.0, 1.0);
    lbl_voigt_lte_matrix_prepare_manylines(mat, flat, atm);

    const int cores = arts_omp_get_max_threads();
    arts_omp_set_num_threads(1);

    lbl_voigt_lte_matrix_prepare_manylines(mat, flat, atm);
    arts_omp_set_num_threads(cores);

    lbl_voigt_lte_matrix_prepare_sort(mat);
  }

  {
    constexpr Index M = 100;
    constexpr Index N = 10'000;
    AbsorptionBands bands{};
    for (Index i = 0; i < N; ++i) {
      QuantumIdentifier key{bnd_qid};
      key.state[QuantumNumberType::v] = {.upper = i, .lower = i};
      bands[key] = AbsorptionBand{.lines = create_lines(M)};
    }
    AtmPoint atm;
    atm.temperature        = 250.0;
    atm.pressure           = 182.0;
    atm[bnd_qid.isot.spec] = 0.21;
    Matrix mat(M, 5);
    const std::vector<lbl::flat_band_data> flat = lbl::flatter_view(
        bands, ZeemanPolarization::no, [](auto&, auto&) { return true; });

    lbl_voigt_lte_matrix_prepare_manybands(mat, flat, atm);

    const int cores = arts_omp_get_max_threads();
    arts_omp_set_num_threads(1);

    lbl_voigt_lte_matrix_prepare_manybands(mat, flat, atm);
    arts_omp_set_num_threads(cores);
  }

  {
    constexpr Index M = 10'000;
    constexpr Index N = 2'000;
    constexpr Index n = 5;
    const Matrix mat  = random_numbers<2>({M, n}, 0.0, 1.0);
    Vector f          = random_numbers<1>({N}, 0.0, 1.0);
    stdr::sort(f);
    ComplexVector a(N);

    buf += lbl_voigt_lte_matrix_sumup(a, mat, f);
    buf += lbl_voigt_lte_matrix_sumup_inf_cutoff(a, mat, f);

    const int cores = arts_omp_get_max_threads();
    arts_omp_set_num_threads(1);

    buf += lbl_voigt_lte_matrix_sumup(a, mat, f);
    buf += lbl_voigt_lte_matrix_sumup_inf_cutoff(a, mat, f);
    arts_omp_set_num_threads(cores);
  }

  {
    constexpr Index M = 1'000'000;
    AbsorptionBands bands{};
    bands[bnd_qid] = AbsorptionBand{.lines = create_lines(M)};
    AtmPoint atm;
    atm.temperature                   = 250.0;
    atm.pressure                      = 182.0;
    atm[bnd_qid.isot.spec]            = 0.21;
    const JacobianTargets jac_targets = create_jac_targets();
    Matrix mat(M, 5 + 5 * jac_targets.target_count());
    const std::vector<lbl::flat_band_data> flat = lbl::flatter_view(
        bands, ZeemanPolarization::no, [](auto&, auto&) { return true; });

    lbl_voigt_lte_matrix_prepare_manylines_jac(mat, flat, atm, jac_targets);

    const int cores = arts_omp_get_max_threads();
    arts_omp_set_num_threads(1);

    lbl_voigt_lte_matrix_prepare_manylines_jac(mat, flat, atm, jac_targets);
    arts_omp_set_num_threads(cores);

    lbl_voigt_lte_matrix_prepare_sort(mat);
  }

  {
    constexpr Index M = 100;
    constexpr Index N = 10'000;
    AbsorptionBands bands{};
    for (Index i = 0; i < N; ++i) {
      QuantumIdentifier key           = bnd_qid;
      key.state[QuantumNumberType::v] = {.upper = i, .lower = i};
      bands[key] = AbsorptionBand{.lines = create_lines(M)};
    }
    AtmPoint atm;
    atm.temperature                   = 250.0;
    atm.pressure                      = 182.0;
    atm[bnd_qid.isot.spec]            = 0.21;
    const JacobianTargets jac_targets = create_jac_targets();
    Matrix mat(M, 5 + 5 * jac_targets.target_count());
    const std::vector<lbl::flat_band_data> flat = lbl::flatter_view(
        bands, ZeemanPolarization::no, [](auto&, auto&) { return true; });

    lbl_voigt_lte_matrix_prepare_manybands_jac(mat, flat, atm, jac_targets);

    const int cores = arts_omp_get_max_threads();
    arts_omp_set_num_threads(1);

    lbl_voigt_lte_matrix_prepare_manybands_jac(mat, flat, atm, jac_targets);
    arts_omp_set_num_threads(cores);
  }

  {
    constexpr Index M = 10'000;
    constexpr Index N = 2'000;
    constexpr Index n = 5 + 5 * 4;
    const Matrix mat  = random_numbers<2>({M, n}, 0.0, 1.0);
    Vector f          = random_numbers<1>({N}, 0.0, 1.0);
    stdr::sort(f);
    ComplexMatrix a(N, n / 5);
    std::vector<bool> df(n / 5 - 1);
    stdr::fill(df, false);

    buf += lbl_voigt_lte_matrix_sumup_jac(a, mat, f, df);
    buf += lbl_voigt_lte_matrix_sumup_inf_cutoff_jac(a, mat, f, df);

    const int cores = arts_omp_get_max_threads();
    arts_omp_set_num_threads(1);

    buf += lbl_voigt_lte_matrix_sumup_jac(a, mat, f, df);
    buf += lbl_voigt_lte_matrix_sumup_inf_cutoff_jac(a, mat, f, df);
    arts_omp_set_num_threads(cores);
  }

  std::println(std::cerr, "Prevent optimizing away: {}", buf);
  arts::print_report();
}