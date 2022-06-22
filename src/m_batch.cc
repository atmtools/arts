/* Copyright (C) 2004-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
                            
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
   USA. */

/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_batch.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2004-09-15 

  \brief  Workspace functions for doing batch calculations.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
using namespace std;

#include "arts.h"
#include "arts_omp.h"
#include "auto_md.h"
#include "math_funcs.h"
#include "physics_funcs.h"
#include "rte.h"
#include "xml_io.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Index GFIELD3_P_GRID;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated 

   2008-07-21 Stefan Buehler */
void ForLoop(Workspace& ws,
             // WS Input:
             const Agenda& forloop_agenda,
             // Control Parameters:
             const Index& start,
             const Index& stop,
             const Index& step,
             const Verbosity& verbosity) {
  CREATE_OUT1;

  for (Index i = start; i <= stop; i += step) {
    out1 << "  Executing for loop body, index: " << i << "\n";
    forloop_agendaExecute(ws, i, forloop_agenda);
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ybatchCalc(Workspace& ws,
                // WS Output:
                ArrayOfVector& ybatch,
                ArrayOfArrayOfVector& ybatch_aux,
                ArrayOfMatrix& ybatch_jacobians,
                // WS Input:
                const Index& ybatch_start,
                const Index& ybatch_n,
                const Agenda& ybatch_calc_agenda,
                // Control Parameters:
                const Index& robust,
                const Verbosity& verbosity) {
  CREATE_OUTS;

  Index first_ybatch_index = 0;

  ArrayOfString fail_msg;
  bool do_abort = false;

  // We allow a start index ybatch_start that is different from 0. We
  // will calculate ybatch_n jobs starting at the start
  // index. Internally, we count from zero, which is the right
  // counting for the output array ybatch. When we call
  // ybatch_calc_agenda, we add ybatch_start to the internal index
  // count.

  // We create a counter, so that we can generate nice output about
  // how many jobs are already done. (All parallel threads later will
  // increment this, so that we really get an accurate total count!)
  Index job_counter = 0;

  // Resize the output arrays:
  ybatch.resize(ybatch_n);
  ybatch_aux.resize(ybatch_n);
  ybatch_jacobians.resize(ybatch_n);
  for (Index i = 0; i < ybatch_n; i++) {
    ybatch[i].resize(0);
    ybatch_aux[i].resize(0);
    ybatch_jacobians[i].resize(0, 0);
  }

  // Go through the batch:

  if (ybatch_n)
#pragma omp parallel for schedule(dynamic) if (!arts_omp_in_parallel() && \
                                               ybatch_n > 1) firstprivate(ws)
    for (Index ybatch_index = first_ybatch_index; ybatch_index < ybatch_n;
         ybatch_index++) {
      Index l_job_counter;  // Thread-local copy of job counter.

      if (do_abort) continue;
#pragma omp critical(ybatchCalc_job_counter)
      { l_job_counter = ++job_counter; }

      {
        ostringstream os;
        os << "  Job " << l_job_counter << " of " << ybatch_n << ", Index "
           << ybatch_start + ybatch_index << ", Thread-Id "
           << arts_omp_get_thread_num() << "\n";
        out2 << os.str();
      }

      try {
        Vector y;
        ArrayOfVector y_aux;
        Matrix jacobian;

        ybatch_calc_agendaExecute(ws,
                                  y,
                                  y_aux,
                                  jacobian,
                                  ybatch_start + ybatch_index,
                                  ybatch_calc_agenda);

        if (y.nelem()) {
#pragma omp critical(ybatchCalc_assign_y)
          ybatch[ybatch_index] = y;
#pragma omp critical(ybatchCalc_assign_y_aux)
          ybatch_aux[ybatch_index] = y_aux;

          // Dimensions of Jacobian:
          const Index Knr = jacobian.nrows();
          const Index Knc = jacobian.ncols();

          if (Knr != 0 || Knc != 0) {
            if (Knr != y.nelem()) {
              ostringstream os;
              os << "First dimension of Jacobian must have same length as the measurement *y*.\n"
                 << "Length of *y*: " << y.nelem() << "\n"
                 << "Dimensions of *jacobian*: (" << Knr << ", " << Knc
                 << ")\n";
              // A mismatch of the Jacobian dimension is a fatal error
              // and should result in program termination. By setting abort
              // to true, this will result in a runtime error in the catch
              // block even if robust == 1
#pragma omp critical(ybatchCalc_setabort)
              do_abort = true;

              throw runtime_error(os.str());
            }

            ybatch_jacobians[ybatch_index] = jacobian;

            // After creation, all individual Jacobi matrices in the array will be
            // empty (size zero). No need for explicit initialization.
          }
        }
      } catch (const std::exception& e) {
        if (robust && !do_abort) {
          ostringstream os;
          os << "WARNING! Job at ybatch_index " << ybatch_start + ybatch_index
             << " failed.\n"
             << "y Vector in output variable ybatch will be empty for this job.\n"
             << "The runtime error produced was:\n"
             << e.what() << "\n";
          out0 << os.str();
        } else {
          // The user wants the batch job to fail if one of the
          // jobs goes wrong.
#pragma omp critical(ybatchCalc_setabort)
          do_abort = true;

          ostringstream os;
          os << "  Job at ybatch_index " << ybatch_start + ybatch_index
             << " failed. Aborting...\n";
          out1 << os.str();
        }
        ostringstream os;
        os << "Run-time error at ybatch_index " << ybatch_start + ybatch_index
           << ": \n"
           << e.what();
#pragma omp critical(ybatchCalc_push_fail_msg)
        fail_msg.push_back(os.str());
      }
    }

  if (fail_msg.nelem()) {
    ostringstream os;

    if (!do_abort) os << "\nError messages from failed batch cases:\n";
    for (ArrayOfString::const_iterator it = fail_msg.begin();
         it != fail_msg.end();
         it++)
      os << *it << '\n';

    if (do_abort)
      throw runtime_error(os.str());
    else
      out0 << os.str();
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ybatchMetProfiles(Workspace& ws,
                       //Output
                       ArrayOfVector& ybatch,
                       //Input
                       const ArrayOfArrayOfSpeciesTag& abs_species,
                       const Agenda& met_profile_calc_agenda,
                       const Vector& f_grid,
                       const Matrix& met_amsu_data,
                       const Matrix& sensor_pos,
                       const Vector& refellipsoid,
                       const Vector& lat_grid,
                       const Vector& lon_grid,
                       const Index& atmosphere_dim,
                       const ArrayOfArrayOfSingleScatteringData& scat_data,
                       //Keyword
                       const Index& nelem_p_grid,
                       const String& met_profile_path,
                       const String& met_profile_pnd_path,
                       const Verbosity& verbosity) {
  GriddedField3 t_field_raw;
  GriddedField3 z_field_raw;
  ArrayOfGriddedField3 vmr_field_raw;
  ArrayOfGriddedField3 pnd_field_raw;
  Vector p_grid;
  Matrix sensor_los;
  Index cloudbox_on;
  ArrayOfIndex cloudbox_limits;
  Matrix z_surface;
  Vector y;
  Index no_profiles = met_amsu_data.nrows();

  //  *vmr_field_raw* is an ArrayOfArrayOfTensor3 where the first array
  //holds the gaseous species.
  //Resize *vmr_field_raw* according to the number of gaseous species
  //elements
  vmr_field_raw.resize(abs_species.nelem());

  //pnd_field_raw is an ArrayOfArrayOfTensor3 where the first array
  //holds the scattering elements.
  // Number of scattering elements:
  const Index N_se = TotalNumberOfElements(scat_data);

  pnd_field_raw.resize(N_se);

  // The satellite zenith angle is read in from the amsu data
  // and converted to arts sensor_los
  ConstVectorView sat_za_from_data = met_amsu_data(Range(joker), 3);

  sensor_los.resize(1, 1);

  // The lat and lon are extracted to get the proper file names of
  // profiles
  ConstVectorView lat = met_amsu_data(Range(joker), 0);
  ConstVectorView lon = met_amsu_data(Range(joker), 1);

  z_surface.resize(1, 1);

  // The spectra .
  y.resize(f_grid.nelem());

  // The batch spectra.
  ybatch.resize(no_profiles);

  // Loop over the number of profiles.
  for (Index i = 0; i < no_profiles; ++i) {
    ostringstream lat_os, lon_os;

    Index lat_prec = 3;
    if (lat[i] < 0) lat_prec--;
    if (abs(lat[i]) >= 10) {
      lat_prec--;
      if (abs(lat[i]) >= 100) lat_prec--;
    }

    lat_os.setf(ios::showpoint | ios::fixed);
    lat_os << setprecision((int)lat_prec) << lat[i];

    Index lon_prec = 4;
    if (lon[i] < 0) lon_prec--;
    if (abs(lon[i]) >= 10) {
      lon_prec--;
      if (abs(lon[i]) >= 100) lon_prec--;
    }
    lon_os.setf(ios::showpoint | ios::fixed);
    lon_os << setprecision((int)lon_prec) << lon[i];

    sensor_los(0, 0) =
        180.0 - (asin(refellipsoid[0] * sin(sat_za_from_data[i] * DEG2RAD) /
                      sensor_pos(0, 0))) *
                    RAD2DEG;

    //Reads the t_field_raw from file
    xml_read_from_file(met_profile_path + "profile.lat_" + lat_os.str() +
                           ".lon_" + lon_os.str() + ".t.xml",
                       t_field_raw,
                       verbosity);

    //Reads the z_field_raw from file
    xml_read_from_file(met_profile_path + "profile.lat_" + lat_os.str() +
                           ".lon_" + lon_os.str() + ".z.xml",
                       z_field_raw,
                       verbosity);

    //Reads the humidity from file - it is only an ArrayofTensor3
    xml_read_from_file(met_profile_path + "profile.lat_" + lat_os.str() +
                           ".lon_" + lon_os.str() + ".H2O.xml",
                       vmr_field_raw[0],
                       verbosity);

    //Reads the pnd_field_raw for one scattering element
    //xml_read_from_file("/rinax/storage/users/rekha/uk_data/profiles/new_obs/newest_forecastfields/reff100/profile.lat_"+lat_os.str()+".lon_"+lon_os.str() + ".pnd100.xml",  pnd_field_raw[0]);

    //xml_read_from_file(met_profile_pnd_path +"reff100_newformat/profile.lat_"+lat_os.str()+".lon_"+lon_os.str() + ".pnd100.xml",  pnd_field_raw[0]);

    xml_read_from_file(met_profile_pnd_path + "lwc_reff15/profile.lat_" +
                           lat_os.str() + ".lon_" + lon_os.str() + ".pnd15.xml",
                       pnd_field_raw[0],
                       verbosity);
    //Write the profile number into a file.
    // xml_write_to_file("profile_number.xml", i);

    // Set z_surface from lowest level of z_field
    z_surface(0, 0) = z_field_raw.data(0, 0, 0);

    /* The vmr_field_raw is an ArrayofArrayofTensor3 where the outer 
     array is for species.
     
     The oxygen and nitrogen VMRs are set to constant values of 0.209
     and 0.782, respectively and are used along with humidity field 
     to generate *vmr_field_raw*.*/

    /*The second element of the species.  The first 3 Tensors in the
    array are the same .  They are pressure grid, latitude grid and
    longitude grid.  The third tensor which is the vmr is set to a 
    constant value of 0.782, corresponding to N2.*/

    vmr_field_raw[1].resize(vmr_field_raw[0]);
    vmr_field_raw[1].copy_grids(vmr_field_raw[0]);
    vmr_field_raw[1] = 0.782;  //vmr of N2

    /*the third element of the species.  the first 3 Tensors in the
    array are the same .  They are pressure grid, latitude grid and
    longitude grid.  The third tensor which is the vmr is set to a 
    constant value of 0.209, corresponding to O2.*/
    vmr_field_raw[2].resize(vmr_field_raw[0]);
    vmr_field_raw[2].copy_grids(vmr_field_raw[0]);
    vmr_field_raw[2] = 0.209;  //vmr of O2

    const Vector& tfr_p_grid =
        t_field_raw.get_numeric_grid(GFIELD3_P_GRID);
    // N_p is the number of elements in the pressure grid
    Index N_p = tfr_p_grid.nelem();

    //Making a p_grid with the first and last element taken from the profile.
    VectorNLogSpace(
        p_grid, nelem_p_grid, tfr_p_grid[0], tfr_p_grid[N_p - 1], verbosity);

    /*To set the cloudbox limits, the lower and upper cloudbox limits
    are to be set.  The lower cloudbox limit is set to the lowest
    pressure level.  The upper level is the highest level where the 
    ice water content is non-zero.*/
    Numeric cl_grid_min, cl_grid_max;

    //Lower limit = lowest pressure level of the original grid.
    //Could it be the interpolated p_grid? FIXME STR
    cl_grid_min = tfr_p_grid[0];

    // A counter for non-zero ice content
    Index level_counter = 0;

    // Loop over all pressure levels
    for (Index ip = 0; ip < N_p; ++ip) {
      //Checking for non-zero ice content. 0.001 is a threshold for
      //ice water content.
      // if((pnd_field_raw[0].data()(ip, 0, 0) > 0.001) || (pnd_field_raw[1](ip, 0, 0) > 0.001))
      if (pnd_field_raw[0].data(ip, 0, 0) > 0.001) {
        ++level_counter;
        //if non-zero ice content is found, it is set to upper
        //cloudbox limit. Moreover, we take one level higher
        // than the upper limit because we want the upper limit
        //to have 0 pnd.
        cl_grid_max = tfr_p_grid[ip + 1];
      }
    }

    //cloudbox limits have dimensions 2*atmosphere_dim
    cloudbox_limits.resize(atmosphere_dim * 2);

    //if there is no cloud in the considered profile, still we
    //need to set the upper limit. I here set the first level
    //for the upper cloudbox limit.
    if (level_counter == 0) {
      cl_grid_max = p_grid[1];
    }

    //Cloudbox is set.
    cloudboxSetManually(cloudbox_on,
                        cloudbox_limits,
                        atmosphere_dim,
                        p_grid,
                        lat_grid,
                        lon_grid,
                        cl_grid_min,
                        cl_grid_max,
                        0,
                        0,
                        0,
                        0,
                        verbosity);

    /*executing the met_profile_calc_agenda
    Agenda communication variables are
    Output of met_profile_calc_agenda : y
    Input to met_profile_calc_agenda  : t_field_raw,
    z_field_raw, vmr_field_raw, pnd_field_raw, p_grid,
    sensor_los, cloudbox_on, cloudbox_limits, z_surface, */

    met_profile_calc_agendaExecute(ws,
                                   y,
                                   t_field_raw,
                                   vmr_field_raw,
                                   z_field_raw,
                                   pnd_field_raw,
                                   p_grid,
                                   sensor_los,
                                   cloudbox_on,
                                   cloudbox_limits,
                                   z_surface,
                                   met_profile_calc_agenda);

    //putting in the spectra *y* for each profile, thus assigning y
    //to the ith row of ybatch
    ybatch[i] = y;

  }  // closing the loop over profile basenames
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ybatchMetProfilesClear(Workspace& ws,
                            //Output
                            ArrayOfVector& ybatch,
                            //Input
                            const ArrayOfArrayOfSpeciesTag& abs_species,
                            const Agenda& met_profile_calc_agenda,
                            const Vector& f_grid,
                            const Matrix& met_amsu_data,
                            const Matrix& sensor_pos,
                            const Vector& refellipsoid,
                            //Keyword
                            const Index& nelem_p_grid,
                            const String& met_profile_path,
                            const Verbosity& verbosity) {
  GriddedField3 t_field_raw;
  GriddedField3 z_field_raw;
  ArrayOfGriddedField3 vmr_field_raw;
  ArrayOfGriddedField3 pnd_field_raw;
  Vector p_grid;
  Matrix sensor_los;
  Index cloudbox_on = 0;
  ArrayOfIndex cloudbox_limits;
  Matrix z_surface;
  Vector y;
  Index no_profiles = met_amsu_data.nrows();
  //Index no_profiles = met_profile_basenames.nelem();
  // The humidity data is stored as  an ArrayOfTensor3 whereas
  // vmr_field_raw is an ArrayOfArrayOfTensor3
  GriddedField3 vmr_field_raw_h2o;

  vmr_field_raw.resize(abs_species.nelem());

  y.resize(f_grid.nelem());
  ybatch.resize(no_profiles);

  Vector sat_za_from_profile;
  sat_za_from_profile = met_amsu_data(Range(joker), 3);
  Numeric sat_za;

  sensor_los.resize(1, 1);

  Vector lat, lon;
  lat = met_amsu_data(Range(joker), 0);
  lon = met_amsu_data(Range(joker), 1);

//   Vector oro_height;
//   oro_height = met_amsu_data(Range(joker), 5);

  z_surface.resize(1, 1);
  for (Index i = 0; i < no_profiles; ++i) {
    ostringstream lat_os, lon_os;

    Index lat_prec = 3;
    if (lat[i] < 0) lat_prec--;
    if (abs(lat[i]) >= 10) {
      lat_prec--;
      if (abs(lat[i]) >= 100) lat_prec--;
    }

    lat_os.setf(ios::showpoint | ios::fixed);
    lat_os << setprecision((int)lat_prec) << lat[i];

    Index lon_prec = 4;
    if (lon[i] < 0) lon_prec--;
    if (abs(lon[i]) >= 10) {
      lon_prec--;
      if (abs(lon[i]) >= 100) lon_prec--;
    }
    lon_os.setf(ios::showpoint | ios::fixed);
    lon_os << setprecision((int)lon_prec) << lon[i];
    cout << lat_os.str() << endl;
    cout << lon_os.str() << endl;

    sat_za = sat_za_from_profile[i];

    sensor_los(Range(joker), 0) =
        180.0 -
        (asin(refellipsoid[0] * sin(sat_za * PI / 180.) / sensor_pos(0, 0))) *
            180. / PI;
    cout << "sensor_los" << sat_za_from_profile[i] << endl;
    cout << "sensor_los" << sat_za << endl;
    cout << "sensor_los" << sensor_los << endl;
    //Reads the t_field_raw from file

    xml_read_from_file(met_profile_path + "profile.lat_" + lat_os.str() +
                           ".lon_" + lon_os.str() + ".t.xml",
                       t_field_raw,
                       verbosity);
    //Reads the z_field_raw from file
    xml_read_from_file(met_profile_path + "profile.lat_" + lat_os.str() +
                           ".lon_" + lon_os.str() + ".z.xml",
                       z_field_raw,
                       verbosity);

    //Reads the humidity from file - it is only an ArrayofTensor3
    // The vmr_field_raw is an ArrayofArrayofTensor3 where the outer
    // array is for species
    xml_read_from_file(met_profile_path + "profile.lat_" + lat_os.str() +
                           ".lon_" + lon_os.str() + ".H2O.xml",
                       vmr_field_raw_h2o,
                       verbosity);
    //xml_read_from_file("/home/home01/rekha/uk/profiles/sat_vmr/profile.lat_"+lat_os.str()//+".lon_"+lon_os.str() + ".H2O_es.xml",
    //                   vmr_field_raw_h2o, verbosity);

    cout
        << "--------------------------------------------------------------------------"
        << endl;
    cout << "The file"
         << met_profile_path + "profile.lat_" + lat_os.str() + ".lon_" +
                lon_os.str()
         << "is executed now" << endl;
    cout
        << "--------------------------------------------------------------------------"
        << endl;
    xml_write_to_file("profile_number.xml", i, FILE_TYPE_ASCII, 0, verbosity);
    // the first element of the species is water vapour.

    // N_p is the number of elements in the pressure grid
    //z_surface(0,0) = oro_height[i]+ 0.01;
    z_surface(0, 0) = z_field_raw.data(0, 0, 0);
    cout << "z_surface" << z_surface << endl;
    const Vector& tfr_p_grid =
        t_field_raw.get_numeric_grid(GFIELD3_P_GRID);
    Index N_p = tfr_p_grid.nelem();

    vmr_field_raw[0] = vmr_field_raw_h2o;

    // the second element of the species.  the first 3 Tensors in the
    //array are the same .  They are pressure grid, latitude grid and
    // longitude grid.  The third tensor which is the vmr is set to a
    // constant value of 0.782.
    vmr_field_raw[1].resize(vmr_field_raw[0]);
    vmr_field_raw[1].copy_grids(vmr_field_raw[0]);
    vmr_field_raw[1].data(joker, joker, joker) = 0.782;

    // the second element of the species.  the first 3 Tensors in the
    //array are the same .  They are pressure grid, latitude grid and
    // longitude grid.  The third tensor which is the vmr is set to a
    // constant value of 0.209.
    vmr_field_raw[2].resize(vmr_field_raw[0]);
    vmr_field_raw[2].copy_grids(vmr_field_raw[0]);
    vmr_field_raw[2].data(joker, joker, joker) = 0.209;

    //xml_write_to_file(met_profile_basenames[i]+ ".N2.xml", vmr_field_raw[1]);
    //xml_write_to_file(met_profile_basenames[i]+ ".O2.xml", vmr_field_raw[2]);

    //Making a p_grid with the first and last element taken from the profile.
    // this is because of the extrapolation problem.

    VectorNLogSpace(
        p_grid, nelem_p_grid, tfr_p_grid[0], tfr_p_grid[N_p - 1], verbosity);
    cout << "t_field_raw[0](0,0,0)" << tfr_p_grid[0] << endl;
    cout << "t_field_raw[0](N_p -1,0,0)" << tfr_p_grid[N_p - 1] << endl;
    xml_write_to_file("p_grid.xml", p_grid, FILE_TYPE_ASCII, 0, verbosity);

    // executing the met_profile_calc_agenda
    met_profile_calc_agendaExecute(ws,
                                   y,
                                   t_field_raw,
                                   vmr_field_raw,
                                   z_field_raw,
                                   pnd_field_raw,
                                   p_grid,
                                   sensor_los,
                                   cloudbox_on,
                                   cloudbox_limits,
                                   z_surface,
                                   met_profile_calc_agenda);

    //putting in the spectra *y* for each profile
    ybatch[i] = y;

  }  // closing the loop over profile basenames
}

/* Workspace method: Doxygen documentation will be auto-generated */
void DOBatchCalc(Workspace& ws,
                 ArrayOfTensor7& dobatch_cloudbox_field,
                 ArrayOfTensor5& dobatch_radiance_field,
                 ArrayOfTensor4& dobatch_irradiance_field,
                 ArrayOfTensor5& dobatch_spectral_irradiance_field,
                 const Index& ybatch_start,
                 const Index& ybatch_n,
                 const Agenda& dobatch_calc_agenda,
                 const Index& robust,
                 const Verbosity& verbosity) {
  CREATE_OUTS;

  Index first_ybatch_index = 0;

  ArrayOfString fail_msg;
  bool do_abort = false;

  // We allow a start index ybatch_start that is different from 0. We
  // will calculate ybatch_n jobs starting at the start
  // index. Internally, we count from zero, which is the right
  // counting for the output array ybatch. When we call
  // ybatch_calc_agenda, we add ybatch_start to the internal index
  // count.

  // We create a counter, so that we can generate nice output about
  // how many jobs are already done. (All parallel threads later will
  // increment this, so that we really get an accurate total count!)
  Index job_counter = 0;

  // Resize the output arrays:
  dobatch_cloudbox_field.resize(ybatch_n);
  dobatch_radiance_field.resize(ybatch_n);
  dobatch_irradiance_field.resize(ybatch_n);
  dobatch_spectral_irradiance_field.resize(ybatch_n);

  // Go through the batch:

  if (ybatch_n)
#pragma omp parallel for schedule(dynamic) if (!arts_omp_in_parallel() && \
                                               ybatch_n > 1) firstprivate(ws)
    for (Index ybatch_index = first_ybatch_index; ybatch_index < ybatch_n;
         ybatch_index++) {
      Index l_job_counter;  // Thread-local copy of job counter.

      if (do_abort) continue;
#pragma omp critical(dobatchCalc_job_counter)
      { l_job_counter = ++job_counter; }

      {
        ostringstream os;
        os << "  Job " << l_job_counter << " of " << ybatch_n << ", Index "
           << ybatch_start + ybatch_index << ", Thread-Id "
           << arts_omp_get_thread_num() << "\n";
        out2 << os.str();
      }

      try {
        Tensor7 cloudbox_field;
        Tensor5 radiance_field;
        Tensor4 irradiance_field;
        Tensor5 spectral_irradiance_field;

        dobatch_calc_agendaExecute(ws,
                                   cloudbox_field,
                                   radiance_field,
                                   irradiance_field,
                                   spectral_irradiance_field,
                                   ybatch_start + ybatch_index,
                                   dobatch_calc_agenda);

#pragma omp critical(dobatchCalc_assign_cloudbox_field)
        dobatch_cloudbox_field[ybatch_index] = cloudbox_field;
#pragma omp critical(dobatchCalc_assign_radiance_field)
        dobatch_radiance_field[ybatch_index] = radiance_field;
#pragma omp critical(dobatchCalc_assign_irradiance_field)
        dobatch_irradiance_field[ybatch_index] = irradiance_field;
#pragma omp critical(dobatchCalc_assign_spectral_irradiance_field)
        dobatch_spectral_irradiance_field[ybatch_index] =
            spectral_irradiance_field;

      } catch (const std::exception& e) {
        if (robust && !do_abort) {
          ostringstream os;
          os << "WARNING! Job at ybatch_index " << ybatch_start + ybatch_index
             << " failed.\n"
             << "element in output variables will be empty for this job.\n"
             << "The runtime error produced was:\n"
             << e.what() << "\n";
          out0 << os.str();
        } else {
          // The user wants the batch job to fail if one of the
          // jobs goes wrong.
#pragma omp critical(dobatchCalc_setabort)
          do_abort = true;

          ostringstream os;
          os << "  Job at ybatch_index " << ybatch_start + ybatch_index
             << " failed. Aborting...\n";
          out1 << os.str();
        }
        ostringstream os;
        os << "Run-time error at ybatch_index " << ybatch_start + ybatch_index
           << ": \n"
           << e.what();
#pragma omp critical(dobatchCalc_push_fail_msg)
        fail_msg.push_back(os.str());
      }
    }

  if (fail_msg.nelem()) {
    ostringstream os;

    if (!do_abort) os << "\nError messages from failed batch cases:\n";
    for (ArrayOfString::const_iterator it = fail_msg.begin();
         it != fail_msg.end();
         it++)
      os << *it << '\n';

    if (do_abort)
      throw runtime_error(os.str());
    else
      out0 << os.str();
  }
}
