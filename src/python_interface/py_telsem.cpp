#include <python_interface.h>

namespace Python {
void py_telsem(py::module_& m) try {
  py::class_<TelsemAtlas>(m, "TelsemAtlas")
      .def("__getstate__", [](TelsemAtlas& self) {
        return py::make_tuple(self.DataCount(),
                              self.ChannelCount(),
                              self.Name(),
                              self.Month(),
                              self.Lat(),
                              self.Cells(),
                              self.FirstCells(),
                              self.Emis(),
                              self.Emis_err(),
                              self.Correlations(),
                              self.Classes1(),
                              self.Classes2(),
                              self.Cellnumber(),
                              self.Correspondance());
      })
      .def("__setstate__", [](TelsemAtlas* self, const py::tuple& t) {
        ARTS_USER_ERROR_IF(t.size() != 14, "Invalid state!")
        
        new (self) TelsemAtlas{};
        self->DataCount() = py::cast<Index>(t[0]);
        self->ChannelCount() = py::cast<Index>(t[1]);
        self->Name() = py::cast<String>(t[2]);
        self->Month() = py::cast<Index>(t[3]);
        self->Lat() = py::cast<Numeric>(t[4]);
        self->Cells() = py::cast<ArrayOfIndex>(t[5]);
        self->FirstCells() = py::cast<ArrayOfIndex>(t[6]);
        self->Emis() = py::cast<Matrix>(t[7]);
        self->Emis_err() = py::cast<Matrix>(t[8]);
        self->Correlations() = py::cast<Tensor3>(t[9]);
        self->Classes1() = py::cast<ArrayOfIndex>(t[10]);
        self->Classes2() = py::cast<ArrayOfIndex>(t[11]);
        self->Cellnumber() = py::cast<ArrayOfIndex>(t[12]);
        self->Correspondance() = py::cast<ArrayOfIndex>(t[13]);
      });
} catch(std::exception& e) {
  throw std::runtime_error(var_string("DEV ERROR:\nCannot initialize telsem\n", e.what()));
}
}  // namespace Python