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
                         const scattering::IrregularZenithAngleGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &);

void xml_read_from_stream(std::istream &is_xml,
                          scattering::IrregularZenithAngleGrid& grid,
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

template <std::floating_point Scalar, scattering::Format fmt, scattering::Representation repr>
void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::PhaseMatrixData<Scalar, fmt, repr>& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &);

template <std::floating_point Scalar, scattering::Format fmt, scattering::Representation repr>
void xml_read_from_stream(std::istream &is_xml,
                          const scattering::PhaseMatrixData<Scalar, fmt, repr>& grid,
                          bifstream *pbifs [[maybe_unused]]);

void xml_write_to_stream(std::ostream &os_xml,
                         const scattering::ZenithAngleGrid& grid,
                         bofstream *pbofs [[maybe_unused]],
                         const String &name);

void xml_read_from_stream(std::istream &is_xml,
                          scattering::ZenithAngleGrid& za_grid,
                          bifstream *pbifs [[maybe_unused]]);


#endif // XML_IO_SCATTERING_H_
