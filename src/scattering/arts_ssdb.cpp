/** \file arts_ssdb.cxx
 *
 * Implementation of ARTS SSDB interface. See arts_ssdb.h.
 *
 * @author Simon Pfreundschuh, 2020
 */
#include "netcdf.h"
#include <scattering/arts_ssdb.h>
#include <filesystem>

namespace scattering {

namespace arts_ssdb {

namespace detail {

////////////////////////////////////////////////////////////////////////////////
// Helper functions.
////////////////////////////////////////////////////////////////////////////////

std::pair<double, double> match_temp_and_freq(std::string group_name) {
  std::regex group_regex(R"(Freq([0-9\.]*)GHz_T([0-9\.]*)K)");
  std::smatch match;
  bool matches = std::regex_match(group_name, match, group_regex);
  if (matches) {
    double freq = std::stod(match[1]);
    double temp = std::stod(match[2]);
    return std::make_pair(freq, temp);
  }
  throw std::runtime_error("Group name doesn't match expected pattern.");
}

std::tuple<bool, double, double, double> match_particle_properties(
    std::filesystem::path path) {
  std::regex file_regex(
      R"(Dmax([0-9]*)um_Dveq([0-9]*)um_Mass([-0-9\.e]*)kg\.nc)");
  std::smatch match;
  std::string filename = path.filename();
  bool matches = std::regex_match(filename, match, file_regex);
  if (matches) {
    double d_max = std::stod(match[1]) * 1e-6;
    double d_eq = std::stod(match[2]) * 1e-6;
    double m = std::stod(match[3]);
    return std::make_tuple(true, d_eq, d_max, m);
  }
  return std::make_tuple(false, 0.0, 0.0, 0.0);
}

std::string match_habit_name(std::filesystem::path path) {
  std::string result = "";
  std::regex folder_regex(R"(([\w]*)_Id([\d]*))");
  std::smatch match;
  for (auto it = path.begin(); it != path.end(); ++it) {
    std::string folder = *it;
    bool matches = std::regex_match(folder, match, folder_regex);
    if (matches) {
      result = match[1];
    }
  }
  return result;
}

void sort_by_d_eq(std::vector<double> &d_eq,
                  std::vector<double> &d_max,
                  std::vector<double> &m) {
  std::vector<size_t> indices{};
  indices.resize(d_eq.size());
  for (size_t i = 0; i < indices.size(); ++i) indices[i] = i;

  auto compare_d_eq = [&d_eq](size_t i, size_t j) { return d_eq[i] < d_eq[j]; };
  std::sort(indices.begin(), indices.end(), compare_d_eq);

  std::vector<double> copy = d_eq;
  for (size_t i = 0; i < indices.size(); ++i) d_eq[i] = copy[indices[i]];
  copy = d_max;
  for (size_t i = 0; i < indices.size(); ++i) d_max[i] = copy[indices[i]];
  copy = m;
  for (size_t i = 0; i < indices.size(); ++i) m[i] = copy[indices[i]];
}
}  // namespace detail

////////////////////////////////////////////////////////////////////////////////
// ScatteringData
////////////////////////////////////////////////////////////////////////////////

void ScatteringData::determine_format() {
  format_ = DataFormat::Gridded;
  if (group_.has_variable("phaMat_data_real")) {
    format_ = DataFormat::Spectral;
  }
  if (group_.has_variable("extMat_data_real")) {
    format_ = DataFormat::FullySpectral;
  }
}

template <typename Float>
math::Vector<Float> ScatteringData::get_vector(std::string name) {
  auto variable = group_.get_variable(name);
  auto size = variable.size();
  auto result = math::Vector<Float>{size};
  variable.read(result.data());
  return result;
}

ParticleType ScatteringData::get_particle_type() {
  if (format_ == DataFormat::Gridded) {
    auto phase_matrix_shape = group_.get_variable("phaMat_data").shape();
    if (phase_matrix_shape[0] == 6) {
      return ParticleType::Random;
    } else {
      return ParticleType::AzimuthallyRandom;
    }
  } else if (format_ == DataFormat::Spectral) {
    auto phase_matrix_shape = group_.get_variable("phaMat_data_real").shape();
    if (phase_matrix_shape[0] == 6) {
      return ParticleType::Random;
    } else {
      return ParticleType::AzimuthallyRandom;
    }
  }
  return ParticleType::AzimuthallyRandom;
}

Index ScatteringData::get_l_max() {
  auto phase_matrix_dimensions =
      group_.get_variable("phaMat_data_real").shape();
  Index l_max = sht::SHT::calc_l_max(phase_matrix_dimensions[3]);
  return l_max;
}

Index ScatteringData::get_n_lon_inc() { return 1; }
Index ScatteringData::get_n_lat_inc() {
  auto variable = group_.get_variable("extMat_data");
  auto dimensions = variable.get_shape_array<math::Index, 3>();
  return dimensions[1];
}
Index ScatteringData::get_n_lon_scat() {
  auto variable = group_.get_variable("phaMat_data");
  auto dimensions = variable.get_shape_array<math::Index, 5>();
  return dimensions[2];
}
Index ScatteringData::get_n_lat_scat() {
  auto variable = group_.get_variable("phaMat_data");
  auto dimensions = variable.get_shape_array<math::Index, 5>();
  return dimensions[4];
}

sht::SHT ScatteringData::get_sht() {
  auto phase_matrix_dimensions =
      group_.get_variable("phaMat_data_real").shape();
  auto l_max = sht::SHT::calc_l_max(phase_matrix_dimensions[3]);
  return sht::SHT(l_max, l_max, 2 * l_max + 2, 2 * l_max + 2);
}

math::Vector<double> ScatteringData::get_lon_inc() {
  if (format_ == DataFormat::Gridded) {
    return get_vector<double>("aa_inc");
  }
  return get_vector<float>("aa_inc").cast<double>();
}

math::Vector<double> ScatteringData::get_lat_inc() {
  if (format_ == DataFormat::Gridded) {
    return get_vector<double>("za_inc");
  }
  return get_vector<float>("za_inc").cast<double>();
}

math::Tensor<double, 7> ScatteringData::get_phase_matrix_data_gridded() {
  // Load data from file.
  auto variable = group_.get_variable("phaMat_data");
  auto dimensions = variable.get_shape_array<math::Index, 5>();
  math::Tensor<double, 5> result{dimensions};
  variable.read(result.data());

  // Reshape and shuffle data.
  math::Tensor<double, 5> result_shuffled = math::cycle_dimensions(result);
  math::Tensor<double, 7> result_reshaped =
      math::unsqueeze<0, 1>(result_shuffled);
  return result_reshaped;
}

math::Tensor<std::complex<double>, 6>
ScatteringData::get_phase_matrix_data_spectral() {
  // Read data from file.
  auto variable_real = group_.get_variable("phaMat_data_real");
  auto variable_imag = group_.get_variable("phaMat_data_imag");
  auto dimensions = variable_real.get_shape_array<math::Index, 4>();
  math::Tensor<float, 4> real{dimensions};
  math::Tensor<float, 4> imag{dimensions};
  variable_real.read(real.data());
  variable_imag.read(imag.data());
  math::Tensor<std::complex<double>, 4> result =
      imag.cast<std::complex<double>>();
  result = result * std::complex<double>(0.0, 1.0);
  result += real;

  // Reshape and shuffle data.
  auto result_shuffled = math::cycle_dimensions(result);
  auto result_reshaped = math::unsqueeze<0, 1>(result_shuffled);
  return result_reshaped;
}

math::Tensor<double, 7> ScatteringData::get_extinction_matrix_data_gridded() {
  // Read data from file.
  auto variable = group_.get_variable("extMat_data");
  auto dimensions = variable.get_shape_array<math::Index, 3>();
  math::Tensor<double, 3> result{dimensions};
  variable.read(result.data());

  // Reshape and shuffle data.
  math::Tensor<double, 3> result_shuffled = math::cycle_dimensions(result);
  math::Tensor<double, 7> result_reshaped =
      math::unsqueeze<0, 1, 4, 5>(result_shuffled);
  return result_reshaped;
}

math::Tensor<std::complex<double>, 6>
ScatteringData::get_extinction_matrix_data_spectral() {
  // Read data from file.
  auto variable = group_.get_variable("extMat_data");
  auto dimensions = variable.get_shape_array<math::Index, 3>();
  math::Tensor<float, 3> result{dimensions};
  variable.read(result.data());

  // Reshape and shuffle data.
  math::Tensor<float, 3> result_shuffled = math::cycle_dimensions(result);
  math::Tensor<float, 6> result_reshaped =
      math::unsqueeze<0, 1, 4>(result_shuffled);
  return result_reshaped.cast<std::complex<double>>();
}

math::Tensor<double, 7> ScatteringData::get_absorption_vector_data_gridded() {
  auto variable = group_.get_variable("absVec_data");
  auto dimensions = variable.get_shape_array<math::Index, 3>();
  math::Tensor<double, 3> result{dimensions};
  variable.read(result.data());

  // Reshape and shuffle data.
  math::Tensor<double, 3> result_shuffled = math::cycle_dimensions(result);
  math::Tensor<double, 7> result_reshaped =
      math::unsqueeze<0, 1, 4, 5>(result_shuffled);
  return result_reshaped;
}

math::Tensor<std::complex<double>, 6>
ScatteringData::get_absorption_vector_data_spectral() {
  auto variable = group_.get_variable("absVec_data");
  auto dimensions = variable.get_shape_array<math::Index, 3>();
  math::Tensor<float, 3> result{dimensions};
  variable.read(result.data());

  // Reshape and shuffle data.
  math::Tensor<float, 3> result_shuffled = math::cycle_dimensions(result);
  math::Tensor<float, 6> result_reshaped =
      math::unsqueeze<0, 1, 4>(result_shuffled);
  return result_reshaped.cast<std::complex<double>>();
}

math::Tensor<double, 7>
ScatteringData::get_backward_scattering_coeff_data_gridded() {
  auto phase_matrix = get_phase_matrix_data_gridded();
  auto dimensions = phase_matrix.dimensions();
  auto backward_scattering_coeff =
      phase_matrix.chip<4>(0).chip<4>(dimensions[5] - 1).chip<4>(0);
  dimensions[4] = 1;
  dimensions[5] = 1;
  dimensions[6] = 1;
  return backward_scattering_coeff.reshape(dimensions);
}

math::Tensor<std::complex<double>, 6>
ScatteringData::get_backward_scattering_coeff_data_spectral() {
  auto phase_matrix = get_phase_matrix_data_spectral();
  auto sht = get_sht();
  auto data_spectral =
      ScatteringDataFieldSpectral(get_f_grid(),
                                  get_t_grid(),
                                  get_lon_inc(),
                                  get_lat_inc(),
                                  sht,
                                  get_phase_matrix_data_spectral());
  auto data_gridded = data_spectral.to_gridded();
  auto phase_matrix_gridded = data_gridded.get_data();
  auto dimensions = phase_matrix_gridded.dimensions();
  auto backward_scattering_coeff =
      phase_matrix_gridded.chip<4>(0).chip<4>(dimensions[5] - 1).chip<4>(0);
  auto dimensions_output = phase_matrix.dimensions();
  dimensions_output[4] = 1;
  dimensions_output[5] = 1;
  return backward_scattering_coeff.cast<std::complex<double>>().reshape(
      dimensions_output);
}

math::Tensor<double, 7>
ScatteringData::get_forward_scattering_coeff_data_gridded() {
  auto phase_matrix = get_phase_matrix_data_gridded();
  auto dimensions = phase_matrix.dimensions();
  auto backward_scattering_coeff =
      phase_matrix.chip<4>(0).chip<4>(0).chip<4>(0);
  dimensions[4] = 1;
  dimensions[5] = 1;
  dimensions[6] = 1;
  return backward_scattering_coeff.reshape(dimensions);
}
math::Tensor<std::complex<double>, 6>
ScatteringData::get_forward_scattering_coeff_data_spectral() {
  auto phase_matrix = get_phase_matrix_data_spectral();
  auto sht = get_sht();

  auto data_spectral =
      ScatteringDataFieldSpectral(get_f_grid(),
                                  get_t_grid(),
                                  get_lon_inc(),
                                  get_lat_inc(),
                                  sht,
                                  get_phase_matrix_data_spectral());
  auto data_gridded = data_spectral.to_gridded();
  auto phase_matrix_gridded = data_gridded.get_data();
  auto forward_scattering_coeff =
      phase_matrix_gridded.chip<4>(0).chip<4>(0).chip<4>(0);
  auto dimensions_output = phase_matrix.dimensions();
  dimensions_output[4] = 1;
  dimensions_output[5] = 1;
  return forward_scattering_coeff.cast<std::complex<double>>().reshape(
      dimensions_output);
}

ScatteringData::operator SingleScatteringDataGridded<double>() {
  assert(format_ == DataFormat::Gridded);

  auto f_grid = std::make_shared<math::Vector<double>>(get_f_grid());
  auto t_grid = std::make_shared<math::Vector<double>>(get_t_grid());
  auto lon_inc = std::make_shared<math::Vector<double>>(get_lon_inc());
  auto lat_inc = std::make_shared<math::Vector<double>>(get_lat_inc());
  auto lon_scat = std::make_shared<math::Vector<double>>(get_lon_scat());
  auto lat_scat = std::make_shared<IrregularLatitudeGrid<double>>(get_lat_scat());
  auto phase_matrix = std::make_shared<math::Tensor<double, 7>>(
      get_phase_matrix_data_gridded());
  auto extinction_matrix = std::make_shared<math::Tensor<double, 7>>(
      get_extinction_matrix_data_gridded());
  auto absorption_vector = std::make_shared<math::Tensor<double, 7>>(
      get_absorption_vector_data_gridded());
  auto backward_scattering_coeff = std::make_shared<math::Tensor<double, 7>>(
      get_backward_scattering_coeff_data_gridded());
  auto forward_scattering_coeff = std::make_shared<math::Tensor<double, 7>>(
      get_backward_scattering_coeff_data_gridded());
  return SingleScatteringDataGridded<double>(f_grid,
                                             t_grid,
                                             lon_inc,
                                             lat_inc,
                                             lon_scat,
                                             lat_scat,
                                             phase_matrix,
                                             extinction_matrix,
                                             absorption_vector,
                                             backward_scattering_coeff,
                                             forward_scattering_coeff);
}

ScatteringData::operator SingleScatteringDataSpectral<double>() {
  assert(format_ == DataFormat::Spectral);

  auto f_grid = std::make_shared<math::Vector<double>>(get_f_grid());
  auto t_grid = std::make_shared<math::Vector<double>>(get_t_grid());
  auto lon_inc = std::make_shared<math::Vector<double>>(get_lon_inc());
  auto lat_inc = std::make_shared<math::Vector<double>>(get_lat_inc());
  auto sht = std::make_shared<sht::SHT>(get_sht());
  auto phase_matrix = std::make_shared<math::Tensor<std::complex<double>, 6>>(
      get_phase_matrix_data_spectral());
  auto extinction_matrix =
      std::make_shared<math::Tensor<std::complex<double>, 6>>(
          get_extinction_matrix_data_spectral());
  auto absorption_vector =
      std::make_shared<math::Tensor<std::complex<double>, 6>>(
          get_absorption_vector_data_spectral());
  auto backward_scattering_coeff =
      std::make_shared<math::Tensor<std::complex<double>, 6>>(
          get_backward_scattering_coeff_data_spectral());
  auto forward_scattering_coeff =
      std::make_shared<math::Tensor<std::complex<double>, 6>>(
          get_backward_scattering_coeff_data_spectral());
  return SingleScatteringDataSpectral<double>(f_grid,
                                              t_grid,
                                              lon_inc,
                                              lat_inc,
                                              sht,
                                              phase_matrix,
                                              extinction_matrix,
                                              absorption_vector,
                                              backward_scattering_coeff,
                                              forward_scattering_coeff);
}

ScatteringData::operator SingleScatteringData() {
  SingleScatteringDataImpl *data = nullptr;
  if (format_ == DataFormat::Gridded) {
    data = new SingleScatteringDataGridded<double>(*this);
  } else if (format_ == DataFormat::Spectral) {
    data = new SingleScatteringDataSpectral<double>(*this);
  }
  return SingleScatteringData(data);
}

////////////////////////////////////////////////////////////////////////////////
// ParticleFile
////////////////////////////////////////////////////////////////////////////////

void ParticleFile::parse_temps_and_freqs() {
  auto group_names = file_handle_.get_group_names();
  std::set<double> freqs;
  std::set<double> temps;
  for (auto &name : group_names) {
    auto freq_and_temp = detail::match_temp_and_freq(name);
    freqs.insert(std::get<0>(freq_and_temp));
    temps.insert(std::get<1>(freq_and_temp));
    auto group = file_handle_.get_group(name).get_group("SingleScatteringData");
    group_map_[freq_and_temp] = group;
  }
  freqs_.resize(freqs.size());
  std::copy(freqs.begin(), freqs.end(), freqs_.begin());
  temps_.resize(temps.size());
  std::copy(temps.begin(), temps.end(), temps_.begin());
  std::sort(freqs_.begin(), freqs_.end());
  std::sort(temps_.begin(), temps_.end());
}

std::array<math::Vector<double>, 4> ParticleFile::get_angular_grids_gridded() {
  std::array<math::Vector<double>, 4> result{};

  Index n_lon_inc_max = 0;
  Index n_lat_inc_max = 0;
  Index n_lon_scat_max = 0;
  Index n_lat_scat_max = 0;

  for (size_t i = 0; i < freqs_.size(); ++i) {
    for (size_t j = 0; j < temps_.size(); ++j) {
      auto data = get_scattering_data(i, j);

      Index n_lon_inc = data.get_n_lon_inc();
      if (n_lon_inc > n_lon_inc_max) {
        result[0] = data.get_lon_inc();
        n_lon_inc_max = n_lon_inc;
      }

      Index n_lat_inc = data.get_n_lat_inc();
      if (n_lat_inc > n_lat_inc_max) {
        result[1] = data.get_lat_inc();
        n_lat_inc_max = n_lat_inc;
      }

      Index n_lon_scat = data.get_n_lon_scat();
      if (n_lon_scat > n_lon_scat_max) {
        result[2] = data.get_lon_scat();
        n_lon_scat_max = n_lon_scat;
      }

      Index n_lat_scat = data.get_n_lat_scat();
      if (n_lat_scat > n_lat_scat_max) {
        result[3] = data.get_lat_scat();
        n_lat_scat_max = n_lat_scat;
      }
    }
  }
  return result;
}

std::tuple<math::Vector<double>, math::Vector<double>, Index>
ParticleFile::get_angular_grids_spectral() {
  std::tuple<math::Vector<double>, math::Vector<double>, Index> result{};

  Index n_lon_inc_max = 0;
  Index n_lat_inc_max = 0;
  Index l_max_max = 0;

  for (size_t i = 0; i < freqs_.size(); ++i) {
    for (size_t j = 0; j < temps_.size(); ++j) {
      auto data = get_scattering_data(i, j);

      Index n_lon_inc = data.get_n_lon_inc();
      if (n_lon_inc > n_lon_inc_max) {
          std::get<0>(result) = data.get_lon_inc();
        n_lon_inc_max = n_lon_inc;
      }

      Index n_lat_inc = data.get_n_lat_inc();
      if (n_lat_inc > n_lat_inc_max) {
          std::get<1>(result) = data.get_lat_inc();
        n_lat_inc_max = n_lat_inc;
      }

      Index l_max = data.get_l_max();
      if (l_max > l_max_max) {
        l_max_max = l_max;
        std::get<2>(result) = l_max;
      }
    }
  }
  return result;
}

ParticleFile::ParticleFile(std::string filename)
    : file_handle_(netcdf4::File::open(filename)) {
  auto properties = detail::match_particle_properties(filename);
  habit_name_ = detail::match_habit_name(filename);
  d_eq_ = std::get<1>(properties);
  d_max_ = std::get<2>(properties);
  mass_ = std::get<3>(properties);
  parse_temps_and_freqs();
}

ParticleType ParticleFile::get_particle_type() {
  auto f = freqs_[0];
  auto t = temps_[0];
  return ScatteringData(group_map_[std::make_pair(f, t)]).get_particle_type();
}

ScatteringData ParticleFile::get_scattering_data(size_t f_index,
                                                 size_t t_index) {
  double freq = freqs_[f_index];
  double temp = temps_[t_index];
  auto found = group_map_.find(std::make_pair(freq, temp));
  return ScatteringData(found->second);
}

ParticleFile::operator SingleScatteringData() {
  auto f_grid = get_f_grid();
  auto t_grid = get_t_grid();

  auto first = ScatteringData(group_map_[std::make_pair(f_grid[0], t_grid[0])]);

  SingleScatteringData result(nullptr);

  if (first.get_format() == DataFormat::Gridded) {
    auto grids = get_angular_grids_gridded();
    auto lon_inc = grids[0] / 180.0 * M_PI;
    auto lat_inc = grids[1] / 180.0 * M_PI;
    auto lon_scat = grids[2] / 180.0 * M_PI;
    auto lat_scat = grids[3] / 180.0 * M_PI;

    result = SingleScatteringData(f_grid,
                                  t_grid,
                                  lon_inc,
                                  lat_inc,
                                  lon_scat,
                                  lat_scat,
                                  first.get_particle_type());
  } else {
    math::Vector<double> lon_inc;
    math::Vector<double> lat_inc;
    Index l_max;

    std::tie(lon_inc, lat_inc, l_max) = get_angular_grids_spectral();
    lon_inc = lon_inc / 180.0 * M_PI;
    lat_inc = lat_inc / 180.0 * M_PI;
    result = SingleScatteringData(f_grid,
                                  t_grid,
                                  lon_inc,
                                  lat_inc,
                                  l_max,
                                  first.get_particle_type());
  }
  for (size_t i = 0; i < freqs_.size(); ++i) {
    for (size_t j = 0; j < temps_.size(); ++j) {
      auto data = get_scattering_data(i, j);
      result.set_data(i, j, data);
    }
  }
  return result;
}

Particle ParticleFile::to_particle() {
  auto properties =
      ParticleProperties{habit_name_, "ARTS SSDB", "", mass_, d_eq_, d_max_, 0.0};
  return Particle(properties, to_single_scattering_data());
}

ParticleFile::DataIterator ParticleFile::begin() {
  return DataIterator(this, 0, 0);
}

ParticleFile::DataIterator ParticleFile::end() {
  return DataIterator(this, freqs_.size(), 0);
}

ParticleFile::DataIterator::DataIterator(const ParticleFile *file,
                                         size_t f_index,
                                         size_t t_index)
    : file_(file), f_index_(f_index), t_index_(t_index) {}

ParticleFile::DataIterator &ParticleFile::DataIterator::operator++() {
  t_index_++;
  if (t_index_ >= file_->temps_.size()) {
    f_index_++;
    t_index_ = 0;
  }
  return *this;
}

ScatteringData ParticleFile::DataIterator::operator*() {
  auto f = file_->freqs_[f_index_];
  auto t = file_->temps_[t_index_];
  return file_->group_map_.find(std::make_pair(f, t))->second;
}

////////////////////////////////////////////////////////////////////////////////
// ParticleFile
////////////////////////////////////////////////////////////////////////////////

void HabitFolder::parse_files() {
  std::vector<double> d_eq_vec, d_max_vec, mass_vec;
  auto it = std::filesystem::directory_iterator(base_path_);
  for (auto &p : it) {
    auto match = detail::match_particle_properties(p.path());
    if (std::get<0>(match)) {
      double d_eq = std::get<1>(match);
      double d_max = std::get<2>(match);
      double mass = std::get<3>(match);
      d_eq_vec.push_back(d_eq);
      d_max_vec.push_back(d_max);
      mass_vec.push_back(mass);
      files_[d_eq] = base_path_ / p.path().filename();
    }
  }
  detail::sort_by_d_eq(d_eq_vec, d_max_vec, mass_vec);
  d_eq_ = math::VectorMap<double>(d_eq_vec.data(), d_eq_vec.size());
  d_max_ = math::VectorMap<double>(d_max_vec.data(), d_max_vec.size());
  mass_ = math::VectorMap<double>(mass_vec.data(), mass_vec.size());
}

HabitFolder::operator ParticleHabit() {
  std::vector<Particle> particles;
  std::cout << "files :: " << files_.size() << std::endl;
  particles.reserve(files_.size());

  ParticleProperties properties{};
  properties.name = "";
  properties.source = "ARTS SSDB";
  properties.refractive_index = "";

  for (Index i = 0; i < d_eq_.size(); ++i) {
    properties.mass = mass_[i];
    properties.d_eq = d_eq_[i];
    properties.d_max = d_max_[i];
    properties.d_aero = 0.0;
    particles.push_back(Particle(properties, ParticleFile(files_[d_eq_[i]])));
    std::cout << "PARTICLE " << i << std::endl;
  }
  return ParticleHabit(particles);
}

HabitFolder::DataIterator HabitFolder::begin() { return DataIterator(this, 0); }

HabitFolder::DataIterator HabitFolder::end() {
  return DataIterator(this, d_eq_.size());
}

}  // namespace arts_ssdb
}  // namespace scattering
