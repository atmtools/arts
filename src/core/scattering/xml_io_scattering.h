#ifndef XML_IO_SCATTERING_H_
#define XML_IO_SCATTERING_H_

#include <xml_io.h>
#include <sht.h>

void xml_read_from_stream(std::istream &is_xml,
                          scattering::sht::SHT &sht,
                          bifstream *pbifs [[maybe_unused]]);

//! Writes SHT to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param sht      SHT
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute (ignored)
 */

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::sht::SHT &sht,
                         bofstream *pbofs [[maybe_unused]],
                         const String &);

void xml_read_from_stream(std::istream &is_xml,
                          scattering::sht::SHT &sht,
                          bifstream *pbifs [[maybe_unused]]);

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::IrregularLatitudeGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &);

void xml_read_from_stream(std::istream &is_xml,
                          scattering::IrregularLatitudeGrid& grid,
                          bifstream *pbifs [[maybe_unused]]);

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::GaussLegendreGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &);

void xml_read_from_stream(std::istream &is_xml,
                          scattering::GaussLegendreGrid& grid,
                          bifstream *pbifs [[maybe_unused]]);

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::DoubleGaussGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &);

void xml_read_from_stream(std::istream &is_xml,
                          scattering::DoubleGaussGrid& grid,
                          bifstream *pbifs [[maybe_unused]]);

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::LobattoGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &);

void xml_read_from_stream(std::istream &is_xml,
                          scattering::LobattoGrid& grid,
                          bifstream *pbifs [[maybe_unused]]);

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::FejerGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &);

void xml_read_from_stream(std::istream &is_xml,
                          scattering::FejerGrid& grid,
                          bifstream *pbifs [[maybe_unused]]);

template <std::floating_point Scalar, scattering::Format fmt, scattering::Representation repr, Index stokes_dim>
void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::PhaseMatrixData<Scalar, fmt, repr, stokes_dim>& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &);

template <std::floating_point Scalar, scattering::Format fmt, scattering::Representation repr, Index stokes_dim>
void xml_read_from_stream(std::istream &is_xml,
                          const scattering::PhaseMatrixData<Scalar, fmt, repr, stokes_dim>& grid,
                          bifstream *pbifs [[maybe_unused]]);

/// Call a function on base class pointer if it can be cast to the given sub-class.
template<typename Base, typename Derived, typename Func>
void call_if_type_matches(Base* base_ptr, Func&& func) {
    if (T* derived_ptr = dynamic_cast<Derived*>(base_ptr)) {
        std::invoke(std::forward<Func>(func), derived_ptr);
    }
}


template<typename Func, typename... Ts>
void callDoSomething(Base* basePtr, Func&& func) {
    (callIfTypeMatches<Ts>(basePtr, std::forward<Func>(func)), ...);
}

#endif // XML_IO_SCATTERING_H_
