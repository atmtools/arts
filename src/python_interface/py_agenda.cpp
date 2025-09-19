#include <nanobind/nanobind.h>
#include <nanobind/stl/bind_map.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>
#include <parameters.h>
#include <workspace.h>
#include <workspace_groups.h>

#include <algorithm>
#include <exception>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "debug.h"
#include "hpy_arts.h"
#include "hpy_vector.h"
#include "python_interface.h"

NB_MAKE_OPAQUE(std::unordered_map<std::string, Wsv>);

extern Parameters parameters;

namespace Python {
Index create_workspace_gin_default_internal(Workspace& ws, const String& key);

std::filesystem::path correct_include_path(
    const std::filesystem::path& path_copy) {
  std::filesystem::path path = path_copy;
  for (auto& prefix : parameters.includepath) {
    if (std::filesystem::is_regular_file(
            path = std::filesystem::path(prefix.c_str()) / path_copy))
      break;
  }

  ARTS_USER_ERROR_IF(not std::filesystem::is_regular_file(path),
                     "Cannot find file: {}\nSearch path(s): {:B,}",
                     path_copy.string(),
                     parameters.includepath)

  // Must add the direcory to include paths as controlfiles know where they are
  std::filesystem::path dir_path = path;
  dir_path.remove_filename();
  String new_inc = dir_path.string();
  if (std::find(parameters.includepath.begin(),
                parameters.includepath.end(),
                new_inc) == parameters.includepath.end())
    parameters.includepath.emplace_back(new_inc);

  return path;
}

void py_agenda(py::module_& m) try {
  py::class_<CallbackOperator> cbd(m, "CallbackOperator");
  generic_interface(cbd);
  cbd.def(
         "__init__",
         [](CallbackOperator* cb,
            const std::function<void(Workspace&)>& f,
            const std::vector<std::string>& i,
            const std::vector<std::string>& o) {
           new (cb)
               CallbackOperator(CallbackOperator::func_t([f](Workspace& ws) {
                                  py::gil_scoped_acquire gil{};
                                  f(ws);
                                }),
                                i,
                                o);
         },
         "f"_a,
         "inputs"_a  = std::vector<std::string>{},
         "outputs"_a = std::vector<std::string>{},
         "Initialize as structured call")
      .def("__call__", [](CallbackOperator& f, Workspace& ws) { f(ws); });

  py::class_<Wsv> wsv(m, "Wsv");
  generic_interface(wsv);
  wsv.def_prop_ro(
         "value",
         [](Wsv& v) { return to_py(v); },
         "A workspace variable.\n\n.. :class:`~pyarts3.arts.Any`")
      .doc() = "A workspace variable wrapper - no manual use required";

  py::class_<Method> methods(m, "Method");
  generic_interface(methods);
  methods
      .def(
          "__init__",
          [](Method* me,
             const std::string& n,
             const std::vector<std::string>& a,
             const std::unordered_map<std::string, std::string>& kw) {
            new (me) Method{n, a, kw};
          },
          "name"_a,
          "args"_a,
          "kwargs"_a,
          "A named method with args and kwargs")
      .def(
          "__init__",
          [](Method* me, const std::string& n, const py::object* const v) {
            new (me) Method{n, from(v).copied()};
          },
          "name"_a,
          "wsv"_a,
          "A method that sets a workspace variable")
      .def_prop_ro(
          "val",
          [](const Method& method) -> py::object {
            const auto& x = method.get_setval();
            if (x) return to_py(x.value());
            return py::none();
          },
          "The value (if any) of a set method.\n\n.. :class:`~pyarts3.arts.Wsv`.\n\n.. :class:`None`")
      .def_prop_ro(
          "name",
          [](const Method& method) { return method.get_name(); },
          "The name of the method.\n\n.. :class:`str`")
      .doc() = "The method class of ARTS";

  auto wsvmap = py::bind_map<std::unordered_map<std::string, Wsv>>(m, "WsvMap");
  wsvmap.doc() = "A helper to speed up file IO of workspace groups";
  wsvmap.def(
      "__init__",
      [](std::unordered_map<std::string, Wsv>* m,
         const std::vector<std::string>& files,
         bool allow_errors) {
        new (m) std::unordered_map<std::string, Wsv>{};
        auto& map = *m;

        std::vector<std::string> actual_error{};
        std::size_t n = internal_workspace_groups().size();

        for (std::size_t i = 0; i < files.size(); ++i) {
          bool happy = false;

          for (std::size_t iwsv = 0; iwsv < n; iwsv++) {
            try {
              auto wsv = Wsv::from_index(iwsv);
              auto f   = wsv.read_from_file(files[i]);
              map.emplace(std::move(f), std::move(wsv));
              happy = true;
              break;
            } catch (std::exception&) {
              // Ignore, try next
            }
          }

          if (not happy) {
            actual_error.push_back(files[i]);
          }
        }

        if (not actual_error.empty()) {
          if (allow_errors) {
            if (map.contains("errors")) return;
            map.emplace("errors", Wsv{actual_error});
          } else {
            throw std::runtime_error(std::format(
                "Could not read {}/{} files as a workspace group:\n{:B,}",
                actual_error.size(),
                files.size(),
                actual_error));
          }
        }
      },
      "files"_a,
      "allow_errors"_a = false,
      R"(Initialize from a list of files

The keys are going to be the file names after the internal
pathing has been applied.

Parameters
----------
files : list of str
    The files to read
allow_errors : bool, optional
    If False (default) then an exception is thrown if any file cannot be read.
    Otherwise the files that cannot be read are put into the workspace variable "errors" as a list of strings.
)");
  wsvmap.def(
      "write_split",
      [](const std::unordered_map<std::string, Wsv>& map,
         std::optional<std::string> basename,
         FileType ftype,
         bool clobber) {
        std::vector<std::string> files;

        for (const auto& [key, value] : map) {
          try {
            String filename = key;
            if (basename) filename = *basename + key;
            value.write_to_file(filename, ftype, clobber);
            files.push_back(std::move(filename));
          } catch (std::exception&) {
            // Ignore write errors
          }
        }

        return files;
      },
      "basename"_a = std::nullopt,
      "ftype"_a    = FileType::ascii,
      "clobber"_a  = true,
      R"(Write the workspace variables to files

The variables are written to the files named by the keys.

Parameters
----------
basename : str, optional
    The base name to use for the output files. The key is appended to this to form the file name.
    If not provided (default) the key is used as the file name.
    Note that it is not ``basename`` + ``"/"`` ``key``, but just ``basename + key``, so if you want
    a directory separator you must provide it.
ftype : FileType, optional
    The file type to write, default is ascii.
clobber : bool, optional
    If True (default) overwrite existing files, if False use a unique file name.

Returns
-------
list of str
    The list of files written.  If basename is provided, the key is appended to it.
    If not, the key is the file name.  If file writing fails an exception is not thrown,
    but the file is not in the returned list.
)");

  py::class_<Agenda> ag(m, "Agenda");
  generic_interface(ag);
  ag.def(py::init<std::string>(), "name"_a, "Create with name")
      .def("document",
           &Agenda::sphinx_list,
           "prep"_a = std::string_view{"- "},
           "Returns a list of methods and state")
      .def("add",
           &Agenda::add,
           "method"_a.none(false),
           R"--(
Adds a method to the Agenda

All workspace variables are defaulted, and all GIN with defaults
create anonymous workspace variables.  All input that are not 
workspace variables are added to the workspace

The input order takes priority over the named argument order,
so Copy(a, out=b) will not even see the b variable.
)--")
      .def(
          "execute",
          [](Agenda& a, Workspace& ws) { a.execute(ws); },
          "ws"_a,
          "Executes the agenda on the provided workspace")
      .def(
          "finalize",
          [](Agenda& a, bool fix) { a.finalize(fix); },
          "fix"_a = false,
          "Finalize the agenda, making it possible to use it in the workspace")
      .def_prop_ro(
          "name",
          [](const Agenda& agenda) { return agenda.get_name(); },
          "The name of the agenda.\n\n.. :class:`str`")
      .def_prop_ro(
          "methods",
          [](const Agenda& agenda) { return agenda.get_methods(); },
          "The methods of the agenda.\n\n.. :class:`list[~pyarts3.arts.Method]`");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize agendas\n{}", e.what()));
}
}  // namespace Python
