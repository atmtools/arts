#pragma once

/** The definition of xml_io_stream_function_declaration
 *
 * Include only in the single source file that defines the XML IO of the function
 * types, followed by an explicit instantiation of xml_io_stream_function_declaration
 * for each of them.  The structs of the functions must be complete there.
 */

#include "xml_io_stream_functional.h"

template <typename F>
void xml_io_stream_function_declaration<F>::write(std::ostream&    os,
                                                  const F&         f,
                                                  bofstream*       pbofs,
                                                  std::string_view name) {
  xml_io_stream_function<F>::write(os, f, pbofs, name);
}

template <typename F>
void xml_io_stream_function_declaration<F>::read(std::istream& is, F& f, bifstream* pbifs) {
  xml_io_stream_function<F>::read(is, f, pbifs);
}
