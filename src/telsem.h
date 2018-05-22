/* Copyright (C) 2017
   Oliver Lemke <olemke@core-dump.info>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA.
*/

/*!
  \file   telsem.h

  \brief  This file contains the definition of the TELSEM atlas format.

  This file implements an interface to the TELSEM model for land surface
  microwave emissivities:

  F. Aires et al, "A Tool to Estimate Land‐Surface Emissivities at
  Microwave frequencies (TELSEM) for use in numerical weather
  prediction," Quarterly Journal of the Royal Meteorological
  Society, vol. 137, (656), pp. 690-699, 2011.

*/

#ifndef telsem_h
#define telsem_h

#include <array>
#include "array.h"
#include "matpackIII.h"
#include "mystring.h"
#include "xml_io.h"


/** A telsem atlas
 *
 * Represents a Telsem2 atlas containing land surface microwave emissivities.
 * Since the Atlas contains emissivities only for land surfaces, the data is
 * stored in a sparse format.
 *
 * The emissivities are represented on an equal area grid and numbered
 * sequentially starting with the first latitude band at -90 degrees and
 * moving up to 90 degrees.
 *
 * The correspondance array contains the data indices for each cellnumber
 * if it is contained in the Atlas and NAN otherwise.
 */
class TelsemAtlas
{

public:

    TelsemAtlas() = default;
    TelsemAtlas(const TelsemAtlas &)  = default;
    TelsemAtlas(      TelsemAtlas &&) = default;
    TelsemAtlas &operator=(const TelsemAtlas &)  = default;
    TelsemAtlas &operator=(      TelsemAtlas &&) = default;
    ~TelsemAtlas() = default;

    /*! Create and read atlas from file.
     *
     * @param filename The path of the file from which to read the atlas
     */
    TelsemAtlas(String filename);

    void  set_month(Index m) {month = m;}
    Index get_month() const {return month;}

    const Tensor3 & get_correl() const {return correl;}
    void set_correl(const Tensor3 &t) {correl = t;}

    /*! Check if cell is contained in atlas.
     *
     * @param cellnumber The cellnumber for given coordinates obtained from
     * calc cellnum.
     */
    bool contains(Index cellnumber) const
    {
        if (cellnumber >= correspondence.nelem()) {
            return false;
        }
        return correspondence[cellnumber] >= 0;
    }

    /*! Class of a given cell.
     * Return the class1 value of the given cell. Indexed by cellnumber
     * obtained from calc_cellnum(...).
     *
     * Throws a runtime error if cellnumber is not contained in the atlas,
     * i.e. is not over land.
     *
     * @param cellnumber The cellnumber for which to lookup the class.
     * @return The index representing the class1 value of the cell.
     */
    Index get_class1(Index cellnumber) const
    {
        Index ind = correspondence[cellnumber];
        if (ind < 0) {
            throw std::runtime_error("The cell is not contained in the atlas.");
        }
        return classes1[ind];
    }

    /*! Class of a given cell.
     * Return the class2 value of the given cell. Indexed by cellnumber
     * obtained from calc_cellnum(...).
     *
     * Throws a runtime error if cellnumber is not contained in the atlas,
     * i.e. is not over land.
     *
     * @param cellnumber The cellnumber for which to lookup the class.
     * @return The index representing the class2 value of the cell.
     */
    Index get_class2(Index cellnumber) const
    {
        Index ind = correspondence[cellnumber];
        if (ind < 0) {
            throw std::runtime_error("The cell is not contained in the atlas.");
        }
        return classes2[ind];
    }

    /*! Verically polarized emissivities at 19, 37 and 85 GHz.
     *
     * The vertically polarized emissivities that are used for the Telsem2
     * emissivity interpolation. The index here is the cellnumber obtained
     * for given latitude and longitude using calc_cellnum(...).
     *
     * @param cellnum The atlas' cellunmber from which to extract the emissivities.
     * @return 3-element vector containing the emissivities.
     */
    Vector get_emis_v(Index i) const {
        Index ind = correspondence[i];
        Vector e_v(3);
        if (ind < 0) {
            throw std::runtime_error("The cell is not contained in the atlas.");
        } else {
            e_v[0] = emis(ind, 0);
            e_v[1] = emis(ind, 3);
            e_v[2] = emis(ind, 5);
        }
        return e_v;
    }

    /*! Horizontally polarized emissivities at 19, 37 and 85 GHz.
     *
     * The horizontally polarized emissivities that are used for the Telsem2
     * emissivity interpolation. The index here is the cellnumber obtained
     * for given latitude and longitude using calc_cellnum(...).
     *
     * @param cellnum The atlas' cellunmber from which to extract the emissivities.
     * @return 3-element vector containing the emissivities.
     */
    Vector get_emis_h(Index cellnum) const {
        Index ind = correspondence[cellnum];
        Vector e_h(3);
        if (ind < 0) {
            throw std::runtime_error("The cell is not contained in the atlas.");
        } else {
            e_h[0] = emis(ind, 1);
            e_h[1] = emis(ind, 4);
            e_h[2] = emis(ind, 6);
        }
        return e_h;
    }

    /*! ConstVectorView on emissivities at given index.
     *
     * The vector containing the seven SSMI emissivities
     * contained in the atlas.
     *
     * Throws a runtime error if cellnumber is not contained in the
     * atlas, i.e. over land.
     *
     * @param The cellnumber
     * @return The ConstVectorView on the emissivities.
     */
    ConstVectorView operator[](Index cellnumber) const {
        Index ind = correspondence[cellnumber];
        if (ind < 0) {
            throw std::runtime_error("The cell is not contained in the atlas.");
        } else {
            return emis(ind, joker);
        }
    }

    /*! Read Telsem Atlas from input stream.
     */
    void read(std::istream &is);

    /*! Compute the number of cells in each latitude band.
    *
    * Telsem surface emissivities are represented on an equal area grid.
    * This functions computes the number of grid cells for each latitude
    * band as well as the cell number of the first grid cell in each
    * band.
    */
    void equare();

    /*! Compute cell indices in sparse data storage.
     *
     * For each cell in the Atlas, this function stores the corresponding
     * index in the data array in the correspondence member of this atlas
     * object.
     */
    void telsem_calc_correspondence();


    /*! Compute the cellnumber corrsponding to given coordinates.
     *
     * Computes the cellnumber of the lon/lat grid that contains the given
     * coordinates.
     *
     * Note: The cell may not be contained in the atlas if
     * it is above the sea. The result should be checked using the
     * contains(...) member function before it is used.
     *
     * @param[in] lat The latitude coordinate for which to compute the
     *                containing cell.
     * @param[out] lon The longitude coordinates for which to compute the
     *                 containing cell.
     */
    Index calc_cellnum(Numeric lat, Numeric lon) const;

    /*! Compute the cellnumber of closest cell contained in atlas.
     *
     * @param[in] lat The latitude coordinate for which to compute the
     *                containing cell.
     * @param[out] lon The longitude coordinates for which to compute the
     *                 containing cell.
     */
    Index calc_cellnum_nearest_neighbor(Numeric lat, Numeric lon) const;

    /*! Compute corrdinates of a given cell.
     *
     * @param[in] cellnum The cell number for which to compute the coordinates.
     * @param[out] The latitiude and longitude coordinates of the cell.
     *
     */
    std::pair<Numeric, Numeric> get_coordinates(Index cellnum) const;

    /*! Inter-/Extrapolate emissivities to given frequency.
     *
     * This function interpolates the SSMI emissivities at 19, 37, 85 GHz to
     * the given frequency value.
     *
     * @param emiss19 The SSMI emissivity at 19 GHz
     * @param emiss37 The SSMI emissivity at 37 GHz
     * @param emiss85 The SSMI emissivity at 85 GHz
     * @param f       The frequency in GHz (!)
     * @param class2  The surface class
     *
     * @return The interpolated emissivity
     */
    Numeric interp_freq2(Numeric emiss19,
                         Numeric emiss37,
                         Numeric emiss85,
                         Numeric f,
                         Index class2) const;

    /*! Interpolat emissivities to given zenith angle and frequency.
     *
     * @param theta The zenith angle
     * @param freq  The frequency in GHz (!!!)
     * @param class1 The surface type class
     * @param class2 The sruface type class
     * @param ev The vertical emissivities from the atlas
     * @param eh The horizontal emissivities from atlas

     * @return A pair containing the interpolated horizontal and vertical
     *         emissivities.
     */
    std::pair<Numeric, Numeric> emis_interp(Numeric theta,
                                            Numeric freq,
                                            Index class1,
                                            Index class2,
                                            const ConstVectorView &ev,
                                            const ConstVectorView &eh) const;

    friend std::ostream& operator<<(std::ostream &os, const TelsemAtlas& ta);
    friend void xml_write_to_stream(ostream &,
                                    const TelsemAtlas&,
                                    bofstream*,
                                    const String&,
                                    const Verbosity&);

    friend void xml_read_from_stream(istream&,
                              TelsemAtlas&,
                              bifstream*,
                              const Verbosity&);
private:

    // Number of lines in the Atlas.
    Index ndat;
    // Number of channels in the Atlas.
    Index nchan;
    // Name of the atlas (including version number).
    String name;
    // Month of the Atlas.
    Index month;
    // Resolution of the Atlas.
    Numeric dlat;
    // Number of cells per lat band.
    ArrayOfIndex ncells;
    // The first cell number of lat band.
    ArrayOfIndex firstcells;
    // Emissivities
    Matrix emis;
    // Emissivity uncertainties.
    Matrix emis_err;
    // Emissivity correlations.
    Tensor3 correl;
    // Surface classes.
    ArrayOfIndex classes1;
    ArrayOfIndex classes2;
    // Cellnumber of each of the pixels in the atlas.
    ArrayOfIndex cellnums;
    // Derived from file data
    ArrayOfIndex correspondence;

    // Regression coefficients.
    static const std::array<Numeric, 30> a0_k0;
    static const std::array<Numeric, 30> a0_k1;
    static const std::array<Numeric, 30> a0_k2;
    static const std::array<Numeric, 30> a0_eveh;
    static const std::array<Numeric, 30> a1_eveh;
    static const std::array<Numeric, 30> a2_eveh;
    static const std::array<Numeric, 30> a3_eveh;
    static const std::array<Numeric, 30> b0_eveh;
    static const std::array<Numeric, 30> b1_eveh;
    static const std::array<Numeric, 30> b2_eveh;
    static const std::array<Numeric, 30> b3_eveh;
    static const std::array<Numeric, 4> rapport43_32;
    static const std::array<Numeric, 4> rapport54_43;
};

typedef Array<TelsemAtlas> ArrayOfTelsemAtlas;

std::ostream& operator<<(std::ostream& os, const TelsemAtlas& ta);

#endif /* telsem_h */

