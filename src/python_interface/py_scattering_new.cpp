#include "py_macros.h"
#include <py_auto_interface.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <scattering/eigen_tensor.h>
#include <scattering/sht.h>
#include <scattering/single_scattering_data.h>
#include <scattering/particle.h>
#include <scattering/particle_habit.h>
#include <scattering/bulk_particle_habit.h>

#include "py_macros.h"

namespace Python {
void py_scattering_new(py::module_& bindings_module) {

  py::module m = bindings_module.def_submodule("scattering_new");

  py::enum_<scattering::DataFormat>(m, "DataFormat")
      .value("Gridded", scattering::DataFormat::Gridded)
      .value("Spectral", scattering::DataFormat::Spectral)
      .value("FullySpectral", scattering::DataFormat::FullySpectral)
      .export_values();

  py::enum_<scattering::ParticleType>(m, "ParticleType")
      .value("Random", scattering::ParticleType::Random)
      .value("AzimuthallyRandom", scattering::ParticleType::AzimuthallyRandom)
      .value("General", scattering::ParticleType::General)
      .export_values();

  ////////////////////////////////////////////////////////////////////////////////
  // Single scattering data
  ////////////////////////////////////////////////////////////////////////////////

  py::class_<scattering::SingleScatteringData>(m, "SingleScatteringData")
      // Constructors
      .def(py::init<>())
      .def(py::init<const scattering::SingleScatteringData &>())
      .def(py::init<scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Tensor<double, 7>,
                    scattering::math::Tensor<double, 7>,
                    scattering::math::Tensor<double, 7>,
                    scattering::math::Tensor<double, 7>,
                    scattering::math::Tensor<double, 7>>())
      .def(py::init<scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::ParticleType>())
      .def(py::init<scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::sht::SHT,
                    scattering::math::Tensor<std::complex<double>, 6>,
                    scattering::math::Tensor<std::complex<double>, 6>,
                    scattering::math::Tensor<std::complex<double>, 6>,
                    scattering::math::Tensor<std::complex<double>, 6>,
                    scattering::math::Tensor<std::complex<double>, 6>>())
      .def(py::init<scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Tensor<std::complex<double>, 6>,
                    scattering::math::Tensor<std::complex<double>, 6>,
                    scattering::math::Tensor<std::complex<double>, 6>,
                    scattering::math::Tensor<std::complex<double>, 6>,
                    scattering::math::Tensor<std::complex<double>, 6>>())
      .def(py::init<scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Vector<double>,
                    scattering::math::Index,
                    scattering::ParticleType>())
      // Class methods
      .def("copy", &scattering::SingleScatteringData::copy)
      .def("get_particle_type",
           &scattering::SingleScatteringData::get_particle_type)
      .def("get_data_format",
           &scattering::SingleScatteringData::get_data_format)
      .def("get_f_grid", &scattering::SingleScatteringData::get_f_grid)
      .def("get_t_grid", &scattering::SingleScatteringData::get_t_grid)
      .def("get_lon_inc", &scattering::SingleScatteringData::get_lon_inc)
      .def("get_lat_inc", &scattering::SingleScatteringData::get_lat_inc)
      .def("get_lon_scat", &scattering::SingleScatteringData::get_lon_scat)
      .def("get_lat_scat", &scattering::SingleScatteringData::get_lat_scat)
      .def("get_n_freqs", &scattering::SingleScatteringData::get_n_freqs)
      .def("get_n_temps", &scattering::SingleScatteringData::get_n_temps)
      .def("get_n_lon_inc", &scattering::SingleScatteringData::get_n_lon_inc)
      .def("get_n_lat_inc", &scattering::SingleScatteringData::get_n_lat_inc)
      .def("get_n_lon_scat", &scattering::SingleScatteringData::get_n_lon_scat)
      .def("get_n_lat_scat", &scattering::SingleScatteringData::get_n_lat_scat)
      .def("get_l_max_scat", &scattering::SingleScatteringData::get_l_max_scat)
      .def("get_m_max_scat", &scattering::SingleScatteringData::get_m_max_scat)
      .def("get_stokes_dim", &scattering::SingleScatteringData::get_stokes_dim)
      .def("set_data", &scattering::SingleScatteringData::set_data)
      .def("interpolate_frequency",
           (scattering::SingleScatteringData(scattering::SingleScatteringData::*)(
               scattering::math::Vector<double>) const) &
           scattering::SingleScatteringData::interpolate_frequency)
      .def("interpolate_temperature",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(
               scattering::math::Vector<double>, bool) const) &
               scattering::SingleScatteringData::interpolate_temperature)
      .def("interpolate_angles",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(
               scattering::math::Vector<double>,
               scattering::math::Vector<double>,
               scattering::math::Vector<double>,
               scattering::math::Vector<double>) const) &
               scattering::SingleScatteringData::interpolate_angles)
      .def("downsample_scattering_angles",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(
               scattering::math::Vector<double>,
               scattering::math::Vector<double>) const) &
               scattering::SingleScatteringData::downsample_scattering_angles)
      .def("get_phase_function",
           &scattering::SingleScatteringData::get_phase_function)
      .def("get_phase_function_spectral",
           &scattering::SingleScatteringData::get_phase_function_spectral)
      .def("get_phase_matrix", (scattering::math::Tensor<double, 8>(scattering::SingleScatteringData::*)(Eigen::Index)const ) &scattering::SingleScatteringData::get_phase_matrix)
      .def("get_phase_matrix", (scattering::math::Matrix<double>(scattering::SingleScatteringData::*)(double, double, double, double, double, double, Eigen::Index)) &scattering::SingleScatteringData::get_phase_matrix)
      .def("get_phase_matrix_data",
           &scattering::SingleScatteringData::get_phase_matrix_data)
      .def("get_phase_matrix_data_spectral",
           &scattering::SingleScatteringData::get_phase_matrix_data_spectral)
      .def("get_extinction_matrix_data",
           &scattering::SingleScatteringData::get_extinction_matrix_data)
      .def("get_extinction_coeff",
           &scattering::SingleScatteringData::get_extinction_coeff)
      .def("get_extinction_matrix", (scattering::math::Tensor<double, 8>(scattering::SingleScatteringData::*)(Eigen::Index)const ) &scattering::SingleScatteringData::get_extinction_matrix)
      .def("get_extinction_matrix", (scattering::math::Matrix<double>(scattering::SingleScatteringData::*)(double, double, double, double, Eigen::Index)const ) &scattering::SingleScatteringData::get_extinction_matrix)
      .def("get_absorption_vector_data",
           &scattering::SingleScatteringData::get_absorption_vector_data)
      .def("get_absorption_coeff",
           &scattering::SingleScatteringData::get_absorption_coeff)
      .def("get_absorption_vector", (scattering::math::Tensor<double, 7>(scattering::SingleScatteringData::*)(Eigen::Index)const ) &scattering::SingleScatteringData::get_absorption_vector)
      .def("get_absorption_vector", (scattering::math::Vector<double>(scattering::SingleScatteringData::*)(double, double, double, double, Eigen::Index)const ) &scattering::SingleScatteringData::get_absorption_vector)
      .def("get_forward_scattering_coeff",
           &scattering::SingleScatteringData::get_forward_scattering_coeff)
      .def("get_backward_scattering_coeff",
           &scattering::SingleScatteringData::get_backward_scattering_coeff)
      .def("__iadd__", &scattering::SingleScatteringData::operator+=)
      .def("__add__", &scattering::SingleScatteringData::operator+)
      .def("__imul__", &scattering::SingleScatteringData::operator*=)
      .def("__mul__", &scattering::SingleScatteringData::operator*)
      .def("normalize", &scattering::SingleScatteringData::normalize)
      .def("regrid", &scattering::SingleScatteringData::regrid)
      .def("set_number_of_scattering_coeffs",
           &scattering::SingleScatteringData::set_number_of_scattering_coeffs)
      .def("to_gridded", &scattering::SingleScatteringData::to_gridded)
      .def("to_spectral",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)() const) &
               scattering::SingleScatteringData::to_spectral)
      .def("to_spectral",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(
               scattering::math::Index, scattering::math::Index) const) &
               scattering::SingleScatteringData::to_spectral)
      .def("to_spectral",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(scattering::math::Index,
                                                    scattering::math::Index,
                                                    scattering::math::Index,
                                                    scattering::math::Index)
                const) &
               scattering::SingleScatteringData::to_spectral)
      .def("to_lab_frame",
           (scattering::SingleScatteringData(
               scattering::SingleScatteringData::*)(scattering::math::Index,
                                                    scattering::math::Index,
                                                    scattering::math::Index)
                const) &
               scattering::SingleScatteringData::to_lab_frame)
      .def("set_stokes_dim", &scattering::SingleScatteringData::set_stokes_dim)
      // Data members
      ;


  ////////////////////////////////////////////////////////////////////////////////
  // Particle
  ////////////////////////////////////////////////////////////////////////////////

  py::class_<scattering::Particle>(m, "Particle")
  // Constructors
  .def_static("from_ssdb", &scattering::Particle::from_ssdb)
  .def(py::init<>())
  .def(py::init<scattering::ParticleProperties,scattering::SingleScatteringData>())
  .def(py::init<double,double,double,scattering::SingleScatteringData>())
  .def(py::init<const scattering::Particle &>())
  // Class methods
  .def("operator=", &scattering::Particle::operator=)
  .def("copy", &scattering::Particle::copy)
  .def("get_name", &scattering::Particle::get_name)
  .def("get_source", &scattering::Particle::get_source)
  .def("get_refractive_index", &scattering::Particle::get_refractive_index)
  .def("get_particle_type", &scattering::Particle::get_particle_type)
  .def("get_data_format", &scattering::Particle::get_data_format)
  .def("get_mass", &scattering::Particle::get_mass)
  .def("get_d_max", &scattering::Particle::get_d_max)
  .def("get_d_eq", &scattering::Particle::get_d_eq)
  .def("get_d_aero", &scattering::Particle::get_d_aero)
  .def("get_f_grid", &scattering::Particle::get_f_grid)
  .def("get_t_grid", &scattering::Particle::get_t_grid)
  .def("get_lon_inc", &scattering::Particle::get_lon_inc)
  .def("get_lat_inc", &scattering::Particle::get_lat_inc)
  .def("get_lon_scat", &scattering::Particle::get_lon_scat)
  .def("get_lat_scat", &scattering::Particle::get_lat_scat)
  .def("interpolate_temperature", (scattering::SingleScatteringData(scattering::Particle::*)(double)const ) &scattering::Particle::interpolate_temperature)
      .def("to_spectral", static_cast<scattering::Particle (scattering::Particle::*)() const>(&scattering::Particle::to_spectral))
  .def("to_lab_frame", (scattering::Particle(scattering::Particle::*)(Eigen::Index, Eigen::Index, Eigen::Index)const ) &scattering::Particle::to_lab_frame)
  .def("regrid", &scattering::Particle::regrid)
  .def("set_stokes_dim", &scattering::Particle::set_stokes_dim)
  .def("needs_t_interpolation", &scattering::Particle::needs_t_interpolation)
  .def("get_data", &scattering::Particle::get_data)
  .def("get_phase_function", &scattering::Particle::get_phase_function)
  .def("get_phase_function_spectral", &scattering::Particle::get_phase_function_spectral)
  .def("get_phase_matrix_data", &scattering::Particle::get_phase_matrix_data)
  .def("get_phase_matrix_data_spectral", &scattering::Particle::get_phase_matrix_data_spectral)
  .def("get_phase_matrix", (scattering::math::Tensor<double, 8>(scattering::Particle::*)(Eigen::Index)const ) &scattering::Particle::get_phase_matrix)
  .def("get_phase_matrix", (scattering::math::Matrix<double>(scattering::Particle::*)(double, double, double, double, double, double, Eigen::Index)) &scattering::Particle::get_phase_matrix)
  .def("get_extinction_coeff", &scattering::Particle::get_extinction_coeff)
  .def("get_extinction_matrix_data", &scattering::Particle::get_extinction_matrix_data)
      .def("get_extinction_matrix", (scattering::math::Tensor<double, 8>(scattering::Particle::*)(Eigen::Index)const ) &scattering::Particle::get_extinction_matrix)
      .def("get_extinction_matrix", (scattering::math::Matrix<double>(scattering::Particle::*)(double, double, double, double, Eigen::Index)) &scattering::Particle::get_extinction_matrix)
  .def("get_absorption_coeff", &scattering::Particle::get_absorption_coeff)
  .def("get_absorption_vector_data", &scattering::Particle::get_absorption_vector_data)
    .def("get_absorption_vector", (scattering::math::Tensor<double, 7>(scattering::Particle::*)(Eigen::Index)const ) &scattering::Particle::get_absorption_vector)
    .def("get_absorption_vector", (scattering::math::Vector<double>(scattering::Particle::*)(double, double, double, double, Eigen::Index)) &scattering::Particle::get_absorption_vector)
  .def("get_forward_scattering_coeff", &scattering::Particle::get_forward_scattering_coeff)
      .def("get_backward_scattering_coeff", &scattering::Particle::get_backward_scattering_coeff)
  .PythonInterfaceBasicRepresentation(scattering::Particle);

  //
  // Particle habit
  //

  py::class_<scattering::ParticleHabit>(m, "ParticleHabit")
      .def_static("liquid_spheres", &scattering::ParticleHabit::liquid_spheres)
      .def("from_ssdb", &scattering::ParticleHabit::from_ssdb)
      .def(py::init([]() {return new scattering::ParticleHabit{}; }))
      .PythonInterfaceBasicRepresentation(scattering::ParticleHabit)
      .def("print", [](const scattering::ParticleHabit& ph) {
          std::cout << "hio " << std::endl;
          std::ostringstream os;
          //os << ph << std::endl;
          return os.str();
      })
      .def("print2", [](const scattering::ParticleHabit& ph) {
          std::cout << "hio " << std::endl;
          std::ostringstream os;
          //os << ph << std::endl;
          return os.str();
      })
  .def(py::pickle(
           [](const scattering::ParticleHabit& v) {
               std::stringstream buffer{};
               v.serialize(buffer);
               return py::make_tuple(py::bytes(buffer.str()));
           },
           [](const py::tuple& t) {
               ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!");
               std::istringstream input{t[0].cast<py::bytes>()};
               return scattering::ParticleHabit::deserialize(input);
           }));

  py::class_<BulkParticleHabit>(m, "BulkParticleHabit")
      .def(py::init([]() {return new BulkParticleHabit{}; }))
      .def(py::init([](const String &name, const ArrayOfSingleScatteringData &scat_data, const ArrayOfScatteringMetaData &meta_data, const Agenda &agenda, const ArrayOfString &pnd_agenda_input){return new BulkParticleHabit{name, scat_data, meta_data, agenda, pnd_agenda_input};}))
      .def(py::pickle(
               [](const BulkParticleHabit& habit) {
                   std::shared_ptr<const Agenda> agenda = habit.get_pnd_agenda();
                   if (agenda.get()) {
                       return py::make_tuple(
                           habit.get_name(),
                           *agenda,
                           habit.get_pnd_agenda_input(),
                           *habit.get_particle_habit(),
                           habit.get_index_start(),
                           habit.get_index_end()
                           );
                   }
                   return py::make_tuple(
                       habit.get_name(),
                       py::none(),
                       habit.get_pnd_agenda_input(),
                       *habit.get_particle_habit(),
                       habit.get_index_start(),
                       habit.get_index_end()
                       );
               },
               [](const py::tuple& t) {
                   ARTS_USER_ERROR_IF(t.size() != 6, "Invalid state!");
                   if (t[1] == py::none()) {
                       return BulkParticleHabit{
                           t[0].cast<std::string>(),
                           t[1].cast<Agenda>(),
                           t[2].cast<ArrayOfString>(),
                           std::make_shared<scattering::ParticleHabit>(t[3].cast<scattering::ParticleHabit>()),
                           };
                   }
                   return BulkParticleHabit{
                       t[0].cast<std::string>(),
                       t[4].cast<Index>(),
                       t[5].cast<Index>(),
                       std::make_shared<scattering::ParticleHabit>(t[3].cast<scattering::ParticleHabit>()),
                       };
                    }));



  py::class_<ScatteringSpecies>(bindings_module, "ScatteringSpecies")
      .def(py::init([]() { return new ScatteringSpecies{}; }))
      .PythonInterfaceCopyValue(ScatteringSpecies)
      .PythonInterfaceWorkspaceVariableConversion(ScatteringSpecies)
      .PythonInterfaceFileIO(ScatteringSpecies)
      .PythonInterfaceBasicRepresentation(ScatteringSpecies)
      .def(py::pickle(
               [](const ScatteringSpecies& species) {
                   auto impl_ptr = species.get();
                   if (impl_ptr) {
                       return py::make_tuple(species.get());
                   } else {
                       return py::make_tuple(py::none());
                   }},
               [](const py::tuple& t) {
                   //ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
                   if (t[0] == py::none()) {
                       return ScatteringSpecies{};
                   }
                   auto habit = std::make_shared<BulkParticleHabit>(t[0].cast<BulkParticleHabit>());
                   return ScatteringSpecies{habit};
               }))
      .PythonInterfaceWorkspaceDocumentation(SingleScatteringData);

  py::class_<Array<ScatteringSpecies>>(bindings_module, "ArrayOfScatteringSpeciesInternal")
      .PythonInterfaceBasicRepresentation(Array<ScatteringSpecies>)
      .PythonInterfaceArrayDefault(ScatteringSpecies);

  py::class_<ArrayOfScatteringSpecies, Array<ScatteringSpecies>>(bindings_module, "ArrayOfScatteringSpecies")
      .PythonInterfaceFileIO(ArrayOfScatteringSpecies)
      .PythonInterfaceCopyValue(ArrayOfScatteringSpecies)
      .PythonInterfaceWorkspaceVariableConversion(ArrayOfScatteringSpecies)
      .PythonInterfaceBasicRepresentation(ArrayOfScatteringSpecies)
      .PythonInterfaceIndexItemAccess(ArrayOfScatteringSpecies)
      .def(py::self == py::self)
      .def(py::init([]() { return new ArrayOfScatteringSpecies{}; }))
      .def(py::init(
          [](Index a, ScatteringSpecies b) { return new ArrayOfScatteringSpecies(a, b); }))
      .def(py::init([](const std::vector<ScatteringSpecies>& v) {
        return new ArrayOfScatteringSpecies{v};
      }))
      .def(
          "append",
          [](ArrayOfScatteringSpecies& x, ScatteringSpecies y) { x.emplace_back(y); },
          py::doc("Appends a ScatteringSpecies at the end of the Array"))
      .def(
          "pop",
          [](ArrayOfScatteringSpecies& x) {
            ScatteringSpecies y = x.back();
            x.pop_back();
            return y;
          },
          py::doc("Pops a ScatteringSpecies from the end of the Array"))
      .def(py::pickle(
          [](const ArrayOfScatteringSpecies& v) {
            auto n = v.size();
            std::vector<ScatteringSpecies> out(n);
            std::copy(v.begin(), v.end(), out.begin());
            return py::make_tuple(std::move(out));
          },
          [](const py::tuple& t) {
            ARTS_USER_ERROR_IF(t.size() != 1, "Invalid state!")
            return new ArrayOfScatteringSpecies{t[0].cast<std::vector<ScatteringSpecies>>()};
          }))
      .PythonInterfaceWorkspaceDocumentation(ArrayOfScatteringSpecies);
  py::implicitly_convertible<std::vector<ScatteringSpecies>, ArrayOfScatteringSpecies>();
}
}  // namespace Python
