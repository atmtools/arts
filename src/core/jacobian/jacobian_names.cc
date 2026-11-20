#include "jacobian_names.h"

#include <xml_io_stream_functional_definition.h>

#include "jac_log.h"
#include "jac_logrel.h"
#include "jac_polyfit.h"
#include "jac_rel.h"
#include "jac_rh.h"

template struct xml_io_stream_function_declaration<Jacobian::atm_vec>;
template struct xml_io_stream_function_declaration<Jacobian::atm_mat>;
template struct xml_io_stream_function_declaration<Jacobian::surf_vec>;
template struct xml_io_stream_function_declaration<Jacobian::surf_mat>;
template struct xml_io_stream_function_declaration<Jacobian::subsurf_vec>;
template struct xml_io_stream_function_declaration<Jacobian::subsurf_mat>;
template struct xml_io_stream_function_declaration<Jacobian::line_vec>;
template struct xml_io_stream_function_declaration<Jacobian::line_mat>;
template struct xml_io_stream_function_declaration<Jacobian::sensor_vec>;
template struct xml_io_stream_function_declaration<Jacobian::sensor_mat>;
template struct xml_io_stream_function_declaration<Jacobian::error_vec>;
template struct xml_io_stream_function_declaration<Jacobian::error_mat>;
