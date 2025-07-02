/** \file
    Header file for work with HITRAN collision induced absorption (CIA). 

    The CIA data are part of the HITRAN distribution. They are described in
    Richard, C., I. E. Gordon, L. S. Rothman, M. Abel, L. Frommhold, M. Gustafsson, 
    J.-M. Hartmann, C. Hermans, W. J. Lafferty, G. S. Orton, K.M. Smith, and H. Tran (2012), 
    New section of the HITRAN database: Collision-induced absorption (CIA),
    J. Quant. Spectrosc. Radiat. Transfer, 113, 1276-1285, doi:10.1016/j.jqsrt.2011.11.004.

    \author Stefan Buehler
    \date   2012-11-30
*/

#ifndef cia_h
#define cia_h

#include <matpack.h>
#include <xml.h>

#include <memory>

#include "mystring.h"
#include "species.h"

// Declare existance of some classes:
class bifstream;
class CIARecord;

using ArrayOfCIARecord = Array<CIARecord>;

std::ostream& operator<<(std::ostream& os, const ArrayOfCIARecord& x);

/* Header with implementation. */
void cia_interpolation(VectorView result,
                       const ConstVectorView& frequency,
                       const Numeric& temperature,
                       const GriddedField2& cia_data,
                       const Numeric& T_extrapolfac,
                       const Index& robust);

Index cia_get_index(const ArrayOfCIARecord& cia_data,
                    const SpeciesEnum sp1,
                    const SpeciesEnum sp2);

CIARecord* cia_get_data(const std::shared_ptr<std::vector<CIARecord>>& cia_data,
                        const SpeciesEnum sp1,
                        const SpeciesEnum sp2);

/** CIA data for a single pair of molecules.
 
 A variable of this class can hold the complete information from one HITRAN CIA file.
 A HITRAN CIA data file can hold several datasets (data for different temperatures but
 fixed frequency range). But all datasets are for the same pair of molecules.
 
 \author Stefan Buehler
 \date   2000-08-21  */
class CIARecord {
 public:
  /** Return each molecule name (as a string) that is associated with this CIARecord.
     
     The CIARecord is defined for a pair of molecules!
     
     \param[in] i Must be either 0 or 1. Then the first or second name of the pair
                  is returned.
     
     */
  [[nodiscard]] String MoleculeName(const Index i) const;

  /** Set each molecule name (from a string) that is associated with this CIARecord.
     
     The CIARecord is defined for a pair of molecules. The molecule names are 
     internally stored as species indices.
     
     \param[in] i Must be either 0 or 1. Then the first or second name of the pair
     is returned.
     \param[in] name The molecule name as a string, e.g., "H2O".
     
     */
  void SetMoleculeName(const Index i, const String& name);

  /** Return CIA species index.
     
     \param[in] i Must be either 0 or 1. Then the first or second species index
     is returned.
     */
  [[nodiscard]] SpeciesEnum Species(const Index i) const;

  /** Return number of datasets in this record.
     */
  [[nodiscard]] Index DatasetCount() const;

  /** Return frequency grid for given dataset.
     */
  [[nodiscard]] ConstVectorView FrequencyGrid(Size dataset) const;

  /** Return temperatur grid for given dataset.
     */
  [[nodiscard]] ConstVectorView TemperatureGrid(Size dataset) const;

  /** Return CIA dataset.
     */
  [[nodiscard]] const GriddedField2& Dataset(Size dataset) const;

  /** Return CIA data.
     */
  [[nodiscard]] const ArrayOfGriddedField2& Data() const;

  /** Return CIA data.
   */
  ArrayOfGriddedField2& Data();

  /** Set CIA species.
     \param[in] first CIA Species.
     \param[in] second CIA Species.
     */
  void SetSpecies(const SpeciesEnum first, const SpeciesEnum second);

  /** Vector version of extract.

     Check whether there is a suitable dataset in the CIARecord and do the 
     interpolation.
     
     \param[out] result CIA value for given frequency grid and temperature.
     \param[in] f_grid Frequency grid.
     \param[in] temperature Scalar temparature.
     \param[in] dataset Index of dataset to use.
     \param[in] robust      Set to 1 to suppress runtime errors (and return NAN values instead).
     */
  void Extract(VectorView result,
               const ConstVectorView& f_grid,
               const Numeric& temperature,
               const Numeric& T_extrapolfac,
               const Index& robust) const;

  /** Scalar version of extract.
     
     Use the vector version, if you can, it is more efficient. This is just a 
     convenience wrapper for it.
     
     \return Scalar CIA value at given frequency and temperature.
     \param[in] frequency Scalar frequency
     \param[in] temperature Scalar temparature
     \param[in] dataset Index of dataset to use 
     \param[in] robust      Set to 1 to suppress runtime errors (and return NAN values instead).
     */
  [[nodiscard]] Numeric Extract(const Numeric& frequency,
                                const Numeric& temperature,
                                const Numeric& T_extrapolfac,
                                const Index& robust) const;

  /** Read CIA catalog file. */
  void ReadFromCIA(const String& filename);

  friend void xml_read_from_stream(std::istream& is_xml,
                                   CIARecord& cr,
                                   bifstream* pbifs);

  /** Append other CIARecord to this. */
  void AppendDataset(const CIARecord& c2);

  [[nodiscard]] std::array<SpeciesEnum, 2> TwoSpecies() const;
  std::array<SpeciesEnum, 2>& TwoSpecies();

  CIARecord() = default;

  CIARecord(ArrayOfGriddedField2 data, SpeciesEnum spec1, SpeciesEnum spec2);

  friend std::ostream& operator<<(std::ostream& os, const CIARecord& cr);

 private:
  /** Append dataset to mdata. */
  void AppendDataset(const Vector& freq,
                     const Vector& temp,
                     const ArrayOfVector& cia);

  /** The data itself, directly from the HITRAN file. 
     
     Dimensions:
     Array dimension: Dataset. One file (one molecule pair) can have
                            different datasets, typically for different temperature
                            or frequency ranges.
     Gridded field dimension 1: Frequency [Hz].
     Gridded field dimension 2: Temperature [K].
     Data: Binary absorption cross-sections in m^5 molec^(-2) 
     
     */
  ArrayOfGriddedField2 mdata;

  /** The pair of molecules associated with these CIA data.
     
     Molecules are specified by their ARTS internal mspecies index! (This has
     to be determined upon reading from a file. Should it ever be written out, it
     has to be mapped to a string again.)
     
     We use a plain C array here, since the length of this is always 2.
     */
  std::array<SpeciesEnum, 2> mspecies;
};

template <>
struct std::formatter<CIARecord> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const CIARecord& v, FmtContext& ctx) const {
    return tags.format(ctx, v.TwoSpecies(), tags.sep(), v.Data());
  }
};

template <>
struct xml_io_stream<CIARecord> {
  static constexpr std::string_view type_name = "CIARecord"sv;

  static void write(std::ostream& os,
                    const CIARecord& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, CIARecord& x, bifstream* pbifs = nullptr);
};

#endif  // cia_h
