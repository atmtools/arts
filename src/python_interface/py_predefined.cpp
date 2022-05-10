#include <predefined/predef.h>
#include <py_auto_interface.h>

#include "predefined/predef_data.h"
#include "py_macros.h"

namespace Python {
void internalMTCKD(
    py::module_& m,
    py::class_<Absorption::PredefinedModel::Hitran::MTCKD::WaterData>& c) {
  c.def(py::init([]() {
     return new Absorption::PredefinedModel::Hitran::MTCKD::WaterData{};
   }))
      .def_readwrite("self_absco_ref",
                     &Absorption::PredefinedModel::Hitran::MTCKD::WaterData::
                         self_absco_ref)
      .def_readwrite(
          "for_absco_ref",
          &Absorption::PredefinedModel::Hitran::MTCKD::WaterData::for_absco_ref)
      .def_readwrite(
          "wavenumbers",
          &Absorption::PredefinedModel::Hitran::MTCKD::WaterData::wavenumbers)
      .def_readwrite(
          "self_texp",
          &Absorption::PredefinedModel::Hitran::MTCKD::WaterData::self_texp);

  m.def(
      "get_foreign_h2o",
      [](const Vector& f,
         Numeric p,
         Numeric t,
         Numeric x,
         PredefinedModelData& data) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::Hitran::MTCKD::compute_foreign_h2o(
            pm, f, p, t, x, data.get<Absorption::PredefinedModel::Hitran::MTCKD::WaterData>());
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::arg("predefined_model_data"),
      py::doc(R"--(Computes foreign absorption using MT CKD Hitran version

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
    predefined_model_data : PredefinedModelData
       As WSV
)--"));

  m.def(
      "get_self_h2o",
      [](const Vector& f,
         Numeric p,
         Numeric t,
         Numeric x,
         PredefinedModelData& data) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::Hitran::MTCKD::compute_self_h2o(
            pm, f, p, t, x, data.get<Absorption::PredefinedModel::Hitran::MTCKD::WaterData>());
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::arg("predefined_model_data"),
      py::doc(R"--(Computes self absorption using MT CKD Hitran version

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
    predefined_model_data : PredefinedModelData
       As WSV
)--"));
}

void internalCKDMT350(py::module_& m) {
  m.def(
      "get_self_h2o",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::CKDMT350::compute_self_h2o(pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::doc(R"--(Computes self absorption using MT CKD version 3.50

Parameters:
    f_grid : :class:`Vector`
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
)--"));

  m.def(
      "get_foreign_h2o",
      [](const Vector& f, Numeric p, Numeric t, Numeric x) -> Vector {
        PropagationMatrix pm(f.nelem());
        Absorption::PredefinedModel::CKDMT350::compute_foreign_h2o(
            pm, f, p, t, x);
        return std::move(pm.Data()).flatten();
      },
      py::arg("f_grid"),
      py::arg("rtp_pressure"),
      py::arg("rtp_temperature"),
      py::arg("x_h2o"),
      py::doc(R"--(Computes foreign absorption using MT CKD version 3.50

Parameters:
    f_grid : Vector
        Frequency grid [Hz]
    rtp_pressure : Numeric
        Pressure value [Pa]
    rtp_temperature : Numeric
        Temperature value [K]
    x_h2o : Numeric
        Ratio of water in the atmosphere in the range [0, 1]
)--"));
}

void py_predefined(py::module_& m) {
  auto predef = m.def_submodule("predef");
  predef.doc() = "Contains predefined absorption models";
  auto CKDMT350 = predef.def_submodule("CKDMT350");
  m.doc() = "Contains implementations for using MT CKD version 3.50";
  auto hitran_mtckd = predef.def_submodule("Hitran").def_submodule("MTCKD");

  py::class_<Absorption::PredefinedModel::Hitran::MTCKD::WaterData>
      hitran_mtckd_data(hitran_mtckd, "WaterData");
  py::class_<PredefinedModelData>(m, "PredefinedModelData")
      .def(py::init([]() { return new PredefinedModelData{}; }))
      .PythonInterfaceWorkspaceVariableConversion(PredefinedModelData)
      .PythonInterfaceFileIO(PredefinedModelData)
      .PythonInterfaceCopyValue(PredefinedModelData)
      .PythonInterfaceBasicRepresentation(PredefinedModelData)
      .def("set",
           [](PredefinedModelData& x,
              Absorption::PredefinedModel::Hitran::MTCKD::WaterData d) {
             x.set(std::move(d));
           })
      .def("get_hitran_mtckd_water_data", [](PredefinedModelData& x) {
        return x.get<Absorption::PredefinedModel::Hitran::MTCKD::WaterData>();
      });

  internalCKDMT350(CKDMT350);
  internalMTCKD(hitran_mtckd, hitran_mtckd_data);
}
}  // namespace Python