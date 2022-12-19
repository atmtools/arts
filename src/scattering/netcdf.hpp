/** An object oriented interface to NetCDF files.
 *
 * This single-include header provides a complete interface
 * to the NetCDF4 c-library.
 *
 * https://github.com/simonpf/netcdfhpp
 *
 * Published under MIT license.
 *
 * Copyright: Simon Pfreundschuh, 2020
 */
#pragma once

#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "netcdf.h"

namespace netcdf4 {
namespace detail {

/** Handle NetCDF error
 *
 * Takes a NetCDF error code returned from a NetCDF library call and
 * throws a std::runtime_error with the error messages if code signals
 * an error.
 *
 * @param error_message Additional error message to display
 *        before the NetCDF error.
 * @param The NetCDF error code returned by the library call.
 */
inline void handle_error(std::string error_message, int error_code) {
  if (error_code != NC_NOERR) {
    std::stringstream error;
    error << error_message << "\n";
    error << nc_strerror(error_code) << std::endl;
    throw std::runtime_error(error.str());
  }
}

/** NetCDF file ID capsule.
 *
 * This wrapper struct manages the lifetime of a netcdf file.
 */
struct FileID {
  ~FileID() { close(); }

  void close() {
    if (open) {
      int error = nc_close(id);
      detail::handle_error("Error closing file: ", error);
      open = false;
    }
  }

  operator int() { return id; }

  int id = 0;
  bool open = false;
};

inline void assert_write_mode(int nc_id) {
  int error = nc_enddef(nc_id);
  if (error != NC_ENOTINDEFINE) {
    detail::handle_error("Error leaving define mode: ", error);
  }
}

inline void assert_define_mode(int nc_id) {
  int error = nc_redef(nc_id);
  if (error != NC_EINDEFINE) {
    detail::handle_error("Error (re)entering define mode: ", error);
  }
}

}  // namespace detail

////////////////////////////////////////////////////////////////////////////////
// Enum constants
////////////////////////////////////////////////////////////////////////////////

enum class CreationMode {
    /// Overwrite existing file.
    Clobber = NC_CLOBBER | NC_NETCDF4,
    /// Avoid overwriting existing file.
    NoClobber = NC_NOCLOBBER | NC_NETCDF4
};

enum class OpenMode {
  Write = NC_WRITE,
  Share = NC_SHARE,
  WriteShare = NC_WRITE | NC_SHARE
};

enum class Type {
  NotAType = NC_NAT,
  Byte = NC_BYTE,
  Char = NC_CHAR,
  Short = NC_SHORT,
  Int = NC_INT,
  Long = NC_LONG,
  Float = NC_FLOAT,
  Double = NC_DOUBLE,
  UByte = NC_UBYTE,
  UShort = NC_USHORT,
  UInt = NC_UINT,
  Int64 = NC_INT64,
  UInt64 = NC_UINT64,
  String = NC_STRING,
  MaxAtomic = NC_MAX_ATOMIC_TYPE
};

inline std::ostream& operator<<(std::ostream& out, const Type& t) {
  std::string name = "";
  switch (t) {
    case Type::NotAType: {
      name = "not_a_type";
      break;
    }
    case Type::Byte: {
      name = "byte";
      break;
    }
    case Type::Char: {
      name = "char";
      break;
    }
    case Type::Short: {
      name = "short";
      break;
    }
    case Type::Int: {
      name = "int";
      break;
    }
    case Type::Float: {
      name = "float";
      break;
    }
    case Type::Double: {
      name = "double";
      break;
    }
    case Type::UByte: {
      name = "unsigned byte";
      break;
    }
    case Type::UShort: {
      name = "unsigned short";
      break;
    }
    case Type::UInt: {
      name = "unsigned int";
      break;
    }
    case Type::Int64: {
      name = "int64";
      break;
    }
    case Type::UInt64: {
      name = "unsigned int64";
      break;
    }
    case Type::String: {
      name = "string";
      break;
    }
  }
  out << name;
  return out;
}

template <typename T>
struct TypeProperties;

template <>
struct TypeProperties<int> {
  static constexpr Type value = Type::Int;
  static constexpr auto write = &nc_put_var_int;
  static constexpr auto write_array = &nc_put_vara_int;
  static constexpr auto write_value = &nc_put_var1_int;
  static constexpr auto write_strided = &nc_put_vars_int;
  static constexpr auto read = &nc_get_var_int;
  static constexpr auto read_array = &nc_get_vara_int;
  static constexpr auto read_value = &nc_get_var1_int;
  static constexpr auto read_strided = &nc_get_vars_int;
};

template <>
struct TypeProperties<float> {
  static constexpr Type value = Type::Float;
  static constexpr auto write = &nc_put_var_float;
  static constexpr auto write_array = &nc_put_vara_float;
  static constexpr auto write_value = &nc_put_var1_float;
  static constexpr auto write_strided = &nc_put_vars_float;
  static constexpr auto read = &nc_get_var_float;
  static constexpr auto read_array = &nc_get_vara_float;
  static constexpr auto read_value = &nc_get_var1_float;
  static constexpr auto read_strided = &nc_get_vars_float;
};

template <>
struct TypeProperties<double> {
  static constexpr Type value = Type::Double;
  static constexpr auto write = &nc_put_var_double;
  static constexpr auto write_array = &nc_put_vara_double;
  static constexpr auto write_value = &nc_put_var1_double;
  static constexpr auto write_strided = &nc_put_vars_double;
  static constexpr auto read = &nc_get_var_double;
  static constexpr auto read_array = &nc_get_vara_double;
  static constexpr auto read_value = &nc_get_var1_double;
  static constexpr auto read_strided = &nc_get_vars_double;
};

template <>
struct TypeProperties<char> {
  static constexpr Type value = Type::Char;
  static constexpr auto write = &nc_put_var_schar;
  static constexpr auto write_array = &nc_put_vara_schar;
  static constexpr auto write_value = &nc_put_var1_schar;
  static constexpr auto write_strided = &nc_put_vars_schar;
  static constexpr auto read = &nc_get_var_schar;
  static constexpr auto read_array = &nc_get_vara_schar;
  static constexpr auto read_value = &nc_get_var1_schar;
  static constexpr auto read_strided = &nc_get_vars_schar;
};

////////////////////////////////////////////////////////////////////////////////
// NetCDF Dimension
////////////////////////////////////////////////////////////////////////////////
/** NetCDF Dimension
 *
 * The Dimension class represents a dimension in a NetCDF file.
 */
struct Dimension {
  Dimension() {}

  /** Create dimension
   *
   * Creates new dimension object setting only its ID. The name
   * and size members can then be used to call nc_inq_dim in order
   * to infer the properties of the dimension.
   *
   * @param id_: The NetCDF ID that identifies the dimension in the
   * file.
   */
  Dimension(int id_) : id(id_) {}

  /** Create new dimension object.
   *
   * Creates new dimension object, that is not yes tied to a dimension
   * in a NetCDF file. This means that setting the id member is still
   * required after the dimension has been successfully created using
   * nc_def_dim.
   *
   * @param name_
   * @param size_
   */
  Dimension(std::string name_, int size_) {
    name_.copy(name, name_.size());
    name[name_.size()] = 0;
    size = size_;
  }

  /// Is this dimension unlimitied?
  bool is_unlimited() { return unlimited; }

  /// The dimension ID that is used by the NetCDF-c library.
  int id = 0;
  /// The size of the dimension.
  size_t size = 0;
  /// Whether or not this is an unlimited dimension.
  bool unlimited = false;
  /// The name of the dimension.
  char name[NC_MAX_NAME + 1] = {0};
};

////////////////////////////////////////////////////////////////////////////////
// NetCDF Variable
////////////////////////////////////////////////////////////////////////////////
/** NetCDF variable
 *
 * Represents a variable in a NetCDF File.
 *
 */
class Variable {

 private:

  // Parses dimensions of this variable.
  void parse_dimensions() {
    auto n_dims = dimensions_.capacity();
    auto dim_ids = std::make_unique<int[]>(n_dims);
    int error = nc_inq_vardimid(parent_id_, id_, dim_ids.get());
    detail::handle_error("Error inquiring dimension IDs of variable:", error);

    for (decltype(n_dims) i = 0; i < n_dims; ++i) {
      int dim_id = dim_ids[i];
      Dimension dim{dim_id};
      error = nc_inq_dim(parent_id_, dim_id, dim.name, &dim.size);
      detail::handle_error("Error inquiring dimension:", error);
      dimensions_.push_back(dim);
    }
  }

  // Checks that NetCDF type is compatible with provided C++ type.
  template <typename T>
  void check_type() {
      if (TypeProperties<T>::value != type_) {
          std::stringstream msg;
          msg << "Provided type " << TypeProperties<T>::value << " is incompatible "
              << "with NetCDF type " << type_ << std::endl;
          throw std::runtime_error(msg.str());
      }
  }

public:

  Variable() {}

  /** Create variable.
   *
   * This creates a variable and parses its dimensions and attributes.
   *
   * @param file_ptr Shared pointer to file object to which this variable belongs.
   * @param parent_id NetCDF ID of the parent group of this variable.
   * @param id The variabl ID that identifies this variable in the NetCDF-c library.
   */
  Variable(std::shared_ptr<detail::FileID> file_ptr, int parent_id, int id)
      : id_(id), parent_id_(parent_id), file_ptr_(file_ptr) {
    int n_dims, n_attrs, type;
    int error = nc_inq_var(parent_id_, id_, name_, &type, &n_dims, 0, &n_attrs);
    detail::handle_error("Error inquiring variable:", error);
    type_ = Type(type);
    dimensions_.reserve(n_dims);
    parse_dimensions();
  }

  /** Write data to variable.
     *
     * @tparam T The datatype to write to the variable.
     * @param data Start pointer to the destination of the read operation.
     */
  template <typename T>
  void write(T* data) {
    using TypeTraits = TypeProperties<T>;
    check_type<T>();
    detail::assert_write_mode(*file_ptr_);
    TypeTraits::write(parent_id_, id_, data);
  }

  /** Write data to variable.
   *
   * Write data to hyperslab of variable memory.
   *
   * @tparam T The datatype to write to the variable.
   * @tparam N_DIM The number of dimensions of the variable.
   * @param starts Array containing the start indices of the hyper-slab specifying
   *     the target for the write operation.
   * @param counts Array containing the lengths of the hyper-slab specifying the
   *     the target for the write operation.
   * @param data Start pointer to the data to write to the variable.
   */
  template <typename T, size_t N_DIMS>
  void write(std::array<size_t, N_DIMS> starts,
             std::array<size_t, N_DIMS> counts,
             const T* data) {
    using TypeTraits = TypeProperties<T>;
    check_type<T>();
    detail::assert_write_mode(*file_ptr_);
    TypeTraits::write_array(
        parent_id_, id_, starts.data(), counts.data(), data);
  }

    /** Write single-valued variable.
     *
     * Write given value to a single-valued variable. If the variable is
     * multi-valued, i.e. a (multi-dimensional) array, the value is
     * written to the first element of the (linearized) array.
     *
     * @tparam The type of the variable.
     * @param The value to write.
     */
    template <typename T>
    void write(T t) {
        using TypeTraits = TypeProperties<T>;
        check_type<T>();
        detail::assert_write_mode(*file_ptr_);
        TypeTraits::write_value(parent_id_, id_, 0, &t);
    }

    /** Read all data from variable.
     *
     * @tparam T The datatype to write to the variable.
     * @param data Start pointer to the destination of the read operation.
     */
    template <typename T>
    void read(T* data) {
        using TypeTraits = TypeProperties<T>;
        check_type<T>();
        detail::assert_write_mode(*file_ptr_);
        TypeTraits::read(parent_id_, id_, data);
    }

  /** Read hyperslab of data from variable.
    *
    * Read data from hyperslab of variable memory.
    *
    * @tparam T The datatype to write to the variable.
    * @tparam N_DIM The number of dimensions of the variable.
    * @param starts Array containing the start indices of the hyper-slab specifying
    *     the source of the read operation.
    * @param counts Array containing the lengths of the hyper-slab specifying the
    *     the source of the read operation.
    * @param data Start pointer to the destination of the read operation.
    */
  template <typename T, size_t N_DIMS>
  void read(std::array<size_t, N_DIMS> starts,
            std::array<size_t, N_DIMS> counts,
            T* data) {
    using TypeTraits = TypeProperties<T>;
    check_type<T>();
    detail::assert_write_mode(*file_ptr_);
    TypeTraits::read_array(parent_id_, id_, starts.data(), counts.data(), data);
  }


  /** Read single-valued variable.
   *
   * Read the value of a single-valued variable and returns it by value. If the
   * variable is multi-valued, i.e. a (multi-dimensional) array, the first element
   * of the (linearized) array is returned.
   *
   * @tparam The type of the variable.
   * @return The value of the variable.
   */
  template <typename T>
  T read() {
      using TypeTraits = TypeProperties<T>;
      check_type<T>();
      T result;
      detail::assert_write_mode(*file_ptr_);
      TypeTraits::read_value(parent_id_, id_, 0, &result);
      return result;
  }

  /// Return reference to dimension vector.
  const std::vector<Dimension>& get_dimensions() const { return dimensions_; }

  /// Total number of elements in the variable's data array.
  size_t size() {
    size_t result = 1;
    for (auto& d : dimensions_) {
      result *= d.size;
    }
    return result;
  }

  /// Array containing the sizes of the variable's data array along each dimension.
  std::vector<size_t> shape() {
    std::vector<size_t> result(dimensions_.size());
    for (size_t i = 0; i < result.size(); ++i) {
      result[i] = dimensions_[i].size;
    }
    return result;
  }

  template <typename Index, size_t N>
  std::array<Index, N> get_shape_array() {
    std::array<Index, N> result = {0};
    for (size_t i = 0; i < N; ++i) {
      result[i] = dimensions_[i].size;
    }
    return result;
  }

  /// The variable's name.
  std::string get_name() const { return name_; }

 private:
  int id_, parent_id_;
  std::vector<Dimension> dimensions_;
  char name_[NC_MAX_NAME + 1] = {0};
  Type type_ = Type::NotAType;
  std::shared_ptr<detail::FileID> file_ptr_ = nullptr;
};

////////////////////////////////////////////////////////////////////////////////
// NetCDF Group
////////////////////////////////////////////////////////////////////////////////
/** A NetCDF group.
  *
  * A NetCDF group contains dimensions, variables, attributes
  * and nested groups.
  *
  */
class Group {
 private:

  // Parses dimensions in this group.
  void parse_dimensions() {
    auto dim_ids = std::make_unique<int[]>(n_dims_);
    int error = nc_inq_dimids(id_, 0, dim_ids.get(), 0);
    detail::handle_error("Error inquiring dimension IDs:", error);
    for (int i = 0; i < n_dims_; ++i) {
      int dim_id = dim_ids[i];
      Dimension dim{dim_id};
      error = nc_inq_dim(id_, dim_id, dim.name, &dim.size);
      detail::handle_error("Error inquiring dimensions", error);
      dimensions_[dim.name] = dim;
    }

    dim_ids = std::make_unique<int[]>(n_unl_dims_);
    error = nc_inq_unlimdims(id_, 0, dim_ids.get());
    detail::handle_error("Error inquiring unlimited dimension IDs:", error);
    for (int i = 0; i < n_unl_dims_; ++i) {
      int dim_id = dim_ids[i];
      Dimension dim{dim_id};
      dim.unlimited = true;
      dim.size = -1;
      error = nc_inq_dim(id_, dim_id, dim.name, &dim.size);
      detail::handle_error("Error inquiring dimensions", error);
      dimensions_[dim.name] = dim;
    }
  }

  // Parses variables in this group.
  void parse_variables() {
    auto var_ids = std::make_unique<int[]>(n_vars_);
    int error = nc_inq_varids(id_, 0, var_ids.get());
    detail::handle_error("Error inquiring variable IDs:", error);
    for (int i = 0; i < n_vars_; ++i) {
      int var_id = var_ids[i];
      Variable var{file_ptr_, id_, var_id};
      variables_[var.get_name()] = var;
    }
  }

    // Parses sub-groups of this group
    void parse_groups() {
        auto group_ids = std::make_unique<int[]>(n_groups_);
        int error = nc_inq_grps(id_, 0, group_ids.get());
        detail::handle_error("Error inquiring group IDs:", error);
        for (int i = 0; i < n_groups_; ++i) {
            int group_id = group_ids[i];
            size_t name_length = 0;
            error = nc_inq_grpname_len(group_id, &name_length);
            detail::handle_error("Error inquiring group name length:", error);
            auto name_bfr = std::make_unique<char[]>(name_length);
            error = nc_inq_grpname(group_id, name_bfr.get());
            detail::handle_error("Error inquiring group name:", error);

            Group group{file_ptr_, group_id, name_bfr.get()};
            groups_[group.get_name()] = group;
        }
    }

  // Ensure that file is in define mode.
  void assert_define_mode() { detail::assert_define_mode(id_); }

  // Ensure that file in write mode.
  void assert_write_mode() { detail::assert_write_mode(id_); }

public:

    Group() {}
  /** Create new group.
    *
    * Creates a new group and parses its dimensions, variables and attributes.
    *
    * @param file_ptr Shared pointer to FileID object managing the NetCDF file.
    * @param id The group ID that identifies the group in the NetCDF-c library.
    * @param name The name of the group.
    */
  Group(std::shared_ptr<detail::FileID> file_ptr, int id, std::string name)
      : file_ptr_(file_ptr), id_(id), name_(name) {
    //
    // Inquire number of dimensions, variables, attributes and groups.
    //
    int error = nc_inq_dimids(id_, &n_dims_, 0, 0);
    detail::handle_error("Error inquiring number of dimensions:", error);
    error = nc_inq_unlimdims(id_, &n_unl_dims_, 0);
    detail::handle_error("Error inquiring number of unlimited dimensions:",
                         error);
    error = nc_inq_varids(id_, &n_vars_, 0);
    detail::handle_error("Error inquiring number of variables:", error);
    error = nc_inq_natts(id_, &n_attrs_);
    detail::handle_error("Error inquiring number of attributes:", error);
    error = nc_inq_grps(id_, &n_groups_, 0);
    detail::handle_error("Error inquiring number of groups:", error);

    parse_dimensions();
    parse_variables();
    parse_groups();
    //parse_attributes();
  }

  /** Add dimension to group.
    *
    * @param name of the dimension.
    * @param size The dimensions size.
    * @return The newly created dimension.
    */
  Dimension add_dimension(std::string name, int size) {
    assert_define_mode();
    Dimension dim{name, size};
    int error = nc_def_dim(id_, name.c_str(), size, &dim.id);
    detail::handle_error("Error creating dimensions: ", error);
    dimensions_[name] = dim;
    sync();
    return dimensions_[name];
  }

  /** Add unlimited dimension to group.
    *
    * @param name of the dimension.
    * @return The newly created dimension.
    */
  Dimension add_dimension(std::string name) {
    assert_define_mode();
    Dimension dim{name, -1};
    dim.unlimited = true;
    int error = nc_def_dim(id_, name.c_str(), NC_UNLIMITED, &dim.id);
    detail::handle_error("Error creating dimensions: ", error);
    sync();
    dimensions_[name] = dim;
    return dimensions_[name];
  }

  /** Add sub-group to group.
     *
     * @param Name of the group to add.
     * @return The newly created group.
     */
  Group add_group(std::string name) {
    assert_define_mode();
    int group_id = 0;
    int error = nc_def_grp(id_, name.c_str(), &group_id);
    detail::handle_error("Error creating group: ", error);
    sync();
    groups_[name] = Group(file_ptr_, group_id, name);
    return groups_[name];
  }

  void sync() {
    assert_write_mode();
    int error = nc_sync(id_);
    detail::handle_error("Error entering define mode: ", error);
  }

  /** Add variable to group
    *
    * @param name of the dimension.
    * @param dimension Vector of dimension names identifying the dimensions
    *    of the variable.
    * @type Type enumer specifying the variable type.
    * @return Variable object representing the newly created variable
    */
  Variable add_variable(std::string name,
                        std::vector<std::string> dimensions,
                        Type type) {
    assert_define_mode();
    auto n_dims = dimensions.size();
    std::vector<int> dim_ids;
    dim_ids.reserve(n_dims);
    for (decltype(n_dims) i = 0; i < n_dims; ++i) {
      auto& d = dimensions[i];
      auto search = dimensions_.find(d);
      if (search == dimensions_.end()) {
        std::stringstream msg;
        msg << "Dimension " << d << " is not defined.";
        throw std::runtime_error(msg.str());
      }
      dim_ids[i] = search->second.id;
    }
    int var_id = 0;
    int error = nc_def_var(id_,
                           name.c_str(),
                           static_cast<int>(type),
                           n_dims,
                           dim_ids.data(),
                           &var_id);
    detail::handle_error("Error defining variable:", error);
    sync();
    variables_[name] = Variable(file_ptr_, id_, var_id);
    return variables_[name];
  }

  /** Retrieve dimension by name.
    *
    * @param name Name of the dimension.
    * @return The dimensions object corresponding to the given name.
    */
  Dimension get_dimension(std::string name) {
    auto found = dimensions_.find(name);
    if (found != dimensions_.end()) {
      return found->second;
    }
    std::stringstream msg;
    msg << "Dimension " << name << " not found in dimensions.";
    throw std::runtime_error(msg.str());
  }

  /** Retrieve variable by name.
    *
    * @param name Name of the variable.
    * @return The variable object corresponding to the given name.
    */
  Variable get_variable(std::string name) {
    auto found = variables_.find(name);
    if (found != variables_.end()) {
      return found->second;
    }
    std::stringstream msg;
    msg << "Variable " << name << " not found in variables.";
    throw std::runtime_error(msg.str());
  }

  /// Check whether group has variable of given name.
  bool has_variable(std::string name) {
    auto found = variables_.find(name);
    return found != variables_.end();
  }

  /// The group name.
  std::string get_name() const { return name_; }

  /** Retrieve group by name.
    *
    * @param name Name of the group.
    * @return The group object corresponding to the given name.
    */
  Group get_group(std::string name) {
      auto found = groups_.find(name);
      if (found != groups_.end()) {
          return found->second;
      }
      std::stringstream msg;
      msg << "Group " << name << " not found in variables.";
      throw std::runtime_error(msg.str());
  }

    /** Get vector of group names.
     * @return Vector containing the names this group's subgroups.
     */
  std::vector<std::string> get_group_names() {
    std::vector<std::string> names{};
    names.reserve(groups_.size());
    for (auto& pair : groups_) {
      names.push_back(pair.second.get_name());
    }
    return names;
  }

  /// Check whether group has subgroup of given name.
  bool has_group(std::string name) {
    auto found = groups_.find(name);
    return found != groups_.end();
  }

 protected:
  std::shared_ptr<detail::FileID> file_ptr_;
  int id_;
  std::string name_;
  int n_dims_, n_vars_, n_groups_, n_attrs_, n_unl_dims_;
  std::map<std::string, Dimension> dimensions_ = {};
  std::map<std::string, Variable> variables_ = {};
  std::map<std::string, Group> groups_ = {};
};

////////////////////////////////////////////////////////////////////////////////
// NetCDF File
////////////////////////////////////////////////////////////////////////////////
/** A NetCDF4 file.
 *
 * The File class represents a NetCDF4 file. It provides
 * an object-oriented interface to the NetCDF-c library.
 *
 */
class File : public Group {
 public:
  /** Create new NetCDF4 file.
     *
     * @param path The file's path
     * @param mode Creation mode defining whether or not to over-
     *        write an existing file.
     * @return File instance representing the newly created file.
     */
  static File create(std::string path,
                     CreationMode mode = CreationMode::Clobber) {
    auto file = std::make_shared<detail::FileID>();
    int error = nc_create(path.c_str(), static_cast<int>(mode), &file->id);
    detail::handle_error("Error creating file: " + path, error);
    file->open = true;
    return File(file);
  }

  /** Open NetCDF4 file.
     *
     * @param path The path to the file to open.
     * @param mode Opening mode defining whether write access
     *        is required.
     * @return File instance representing the opened file.
     */
  static File open(std::string path, OpenMode mode = OpenMode::Write) {
    auto file = std::make_shared<detail::FileID>();
    int error = nc_open(path.c_str(), static_cast<int>(mode), &file->id);
    detail::handle_error("Error opening file: " + path, error);
    file->open = true;
    return File(file);
  }

  File(std::shared_ptr<detail::FileID> file_ptr)
      : Group(file_ptr, *file_ptr, "") {}

  /// Close the file.
  void close() { file_ptr_->close(); }

 private:
};

}  // namespace netcdf4
