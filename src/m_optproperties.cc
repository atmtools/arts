/*!
  \file   m_optproperties.cc
  \author Sreerekha T.R. <rekha@uni-bremen.de>, 
          Claudia Emde <claudia@sat.physik.uni-bremen.de
  \date   Mon Jun 10 11:19:11 2002 
  \brief  This filecontains workspace methods for calculating the optical 
  properties for the radiative transfer. 

  Optical properties are the extinction matrix, absorption vector and
  scattering vector.  The optical properties for the gases can be
  calculated with or without Zeeman effect.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include "arts.h"
#include "exceptions.h"
#include "array.h"
#include "matpackIII.h"
#include "matpackVII.h"
#include "scatproperties.h"
#include "logic.h"
#include "interpolation.h"
#include "messages.h"
#include "xml_io.h"
#include "optproperties.h"
#include "math_funcs.h"
#include "sorting.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;

#define part_type scat_data_raw[i_pt].ptype
#define f_datagrid scat_data_raw[i_pt].f_grid
#define za_datagrid scat_data_raw[i_pt].za_grid
#define aa_datagrid scat_data_raw[i_pt].aa_grid
#define pha_mat_data_raw scat_data_raw[i_pt].pha_mat_data  //CPD: changed from pha_mat_data
#define ext_mat_data_raw scat_data_raw[i_pt].ext_mat_data  //which wouldn't let me play with
#define abs_vec_data_raw scat_data_raw[i_pt].abs_vec_data  //scat_data_mono.


//! Calculate single particle extinction. 
/*!
  This method computes the extinction matrix for a single particle
  type from the amplitude matrix.

  \param ext_mat_spt Output and Input : extinction matrix for a single 
                     particle type.
  \param amp_mat Input : amplitude matrix for each particle type
  \param scat_za_index  Input : local zenith angle
  \param scat_aa_index  Input : local azimuth angle
  \param f_index  Input : frequency index
  \param f_grid  Input : frequency grid
 
*/
void ext_mat_sptCalc(
                     Tensor3& ext_mat_spt,
                     const Tensor6& amp_mat,
                     const Index& scat_za_index,
                     const Index& scat_aa_index,
                     const Index& f_index,
                     const Vector& f_grid)
                     
{
  Index npt = ext_mat_spt.npages();
    
  Index stokes_dim = ext_mat_spt.nrows();
  
  if ((stokes_dim > 4) || (stokes_dim <1)){
    throw runtime_error("The dimension of stokes vector "
                        "can be only 1,2,3, or 4");
  }
  
  if ( ext_mat_spt.ncols() != stokes_dim){
 
    throw runtime_error(" The columns of ext_mat_spt should "
                        "agree to stokes_dim");
  }
    
  if (amp_mat.ncols() != 8)
    throw runtime_error(
                        "Amplitude matrix must have 8 columns.");
  
 
 
  //find out frequency
  Numeric freq = f_grid[f_index];

  for (Index i = 0; i < npt; ++i)
    {
      ConstVectorView amp_coeffs = amp_mat(i,
                                      scat_za_index,
                                      scat_aa_index,
                                      scat_za_index,
                                      scat_aa_index,
                                      Range(joker)
                                      );
      amp2ext(ext_mat_spt(i,Range(joker),Range(joker)),
              amp_coeffs,
              freq);
    }
}

//! Calculate extinction matrix (spt) for the convergence test.  

/*
  This method computes the extinction matrix for a single particle type
  from the amplitude matrix. This extinction matrix contains only 
  extinction due to scattering, not the extinction due to absorption. 
  The function can be used only for testing.

\param ext_mat_spt Output and Input: extinction matrix for a
single particle type.
\param pha_mat_spt Input : phase matrix for single particle type
\param f_grid  Input : frequency grid
\param scat_za_grid  Input : zenith angle grid
\param scat_aa_grid  Input : azimuth angle grid
\param stokes_dim  Input : stokes dimension

*/
void ext_mat_sptScat(
                     Tensor3& ext_mat_spt,
                     const Tensor5& pha_mat_spt,
                     const Vector& scat_za_grid,
                     const Vector& scat_aa_grid
                     )
                    
{

#ifndef NDEBUG
  Index nza = pha_mat_spt.nbooks();
  Index naa = pha_mat_spt.npages(); 
#endif
  Index npt = pha_mat_spt.nshelves();
  Index stokes_dim = pha_mat_spt.nrows();


  assert (is_size(scat_za_grid, nza));
  assert (is_size(scat_aa_grid, naa));

  if ((stokes_dim > 4) || (stokes_dim <1)){
    throw runtime_error("The dimension of stokes vector "
                        "can be only 1,2,3, or 4");
  }
  
  ext_mat_spt.resize(npt, stokes_dim, stokes_dim);
 
  for (Index i = 0; i < npt; ++i)
    {
      ConstTensor4View pha=pha_mat_spt(i,
                                       Range(joker),
                                       Range(joker),
                                       Range(joker),
                                       Range(joker));
      
      // ConstVectorView za_grid = scat_za_grid[scat_za_index];
      //ConstVectorView aa_grid = scat_aa_grid[scat_aa_index];
      
      amp2ext_scat(ext_mat_spt(i,Range(joker), Range(joker)),
                   pha,
                   scat_za_grid,
                   scat_aa_grid);
     
    }      

}



//! Calculates  phase matrix for a single particle type 

/*! 
  The elements of the phase matrix is calculated from the
  elements of the amplitude matrix.

 \param pha_mat_spt  Output and Input: phase matrix for a single particle type
 \param amp_mat  Input : amplitude matrix
 \param scat_za_index  Input : zenith angle index
 \param scat_aa_index  Input : azimuth angle index
   
*/

void pha_mat_sptCalc(
                     Tensor5& pha_mat_spt,
                     const Tensor6& amp_mat,
                     const Index& scat_za_index,
                     const Index& scat_aa_index)
  
{
  Index npt = pha_mat_spt.nshelves();
  Index stokes_dim = pha_mat_spt.nrows();
  
  if ((stokes_dim > 4) || (stokes_dim <1)){
    throw runtime_error("The dimension of stokes vector" 
			"can be only 1,2,3, or 4");
  }
  
  if (pha_mat_spt.ncols() != stokes_dim ){
    throw runtime_error("The tensor pha_mat_spt should"
                        "have 4 columns");
  }
  if (amp_mat.ncols() != 8){
      throw runtime_error("Amplitude matrix must have 8 columns.");
  }
  
  
  for (Index i = 0; i < npt; ++i)
    {
      ConstTensor3View amp_coeffs = amp_mat(i,
					    scat_za_index,
					    scat_aa_index,
					    Range(joker),
					    Range(joker),
					    Range(joker)
					    );
      amp2pha(pha_mat_spt(i, Range(joker), Range(joker), 
                          Range(joker), Range(joker)),
              amp_coeffs); 
    }
  //  xml_write_to_file("pha_mat_spt.xml", pha_mat_spt); 
}



//! Calculates phase matrix for the single particle types.
/*! 
  In this function the phase matrix calculated
  in the laboratory frame. This function is used in the calculation
  of the scattering integral (*scat_fieldCalc*).
  
  The interpolation of the data on the actual frequency is the first 
  step in this function. 
  
  Then the transformation from the database coordinate system to to 
  laboratory coordinate system is done.

  WS Output:
  \param pha_mat_spt Extinction matrix for a single particle type.
  WS Input:
  \param scat_data_raw Raw data containing optical properties.
  \param scat_za_grid Zenith angle grid for scattering calculation.
  \param scat_aa_grid Azimuth angle grid for scattering calculation.
  \param scat_za_index Index specifying the propagation direction.
  \param scat_aa_index  Index specifying the azimuth angle
  \param f_index Index specifying the frequency.
  \param f_grid Frequency grid.

  \author Claudia Emde
  \date   2003-05-14
*/

void pha_mat_sptFromData( // Output:
                         Tensor5& pha_mat_spt,
                         // Input:
                         const ArrayOfSingleScatteringData& scat_data_raw,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid,
                         const Index& scat_za_index, // propagation directions
                         const Index& scat_aa_index,
                         const Index& f_index,
                         const Vector& f_grid,
                         const Tensor4& scat_theta,
                         const ArrayOfArrayOfArrayOfArrayOfGridPos&
                         scat_theta_gps,
                         const Tensor5& scat_theta_itws
                         )
{

 
  out3 << "Calculate *pha_mat_spt* from database\n";

  const Index N_pt = scat_data_raw.nelem();
  const Index stokes_dim = pha_mat_spt.ncols();

  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                         "must be 1,2,3 or 4");
  }
  
  assert( pha_mat_spt.nshelves() == N_pt );
  
  // Phase matrix in laboratory coordinate system. Dimensions:
  // [frequency, za_inc, aa_inc, stokes_dim, stokes_dim]
  Tensor5 pha_mat_data_int;
  

  // Loop over the included particle_types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      // First we have to transform the data from the coordinate system 
      // used in the database (depending on the kind of particle type 
      // specified by *ptype*) to the laboratory coordinate sytem. 
      
      // Frequency interpolation:
     
      // The data is interpolated on one frequency. 
      pha_mat_data_int.resize(pha_mat_data_raw.nshelves(), pha_mat_data_raw.nbooks(),
                              pha_mat_data_raw.npages(), pha_mat_data_raw.nrows(), 
                              pha_mat_data_raw.ncols());

      
      // Gridpositions:
      GridPos freq_gp;
      gridpos(freq_gp, f_datagrid, f_grid[f_index]); 

      // Interpolationweights:
      Vector itw(2);
      interpweights(itw, freq_gp);
     
      for (Index i_za_sca = 0; i_za_sca < pha_mat_data_raw.nshelves() ; i_za_sca++)
        {
          for (Index i_aa_sca = 0; i_aa_sca < pha_mat_data_raw.nbooks(); i_aa_sca++)
            {
              for (Index i_za_inc = 0; i_za_inc < pha_mat_data_raw.npages(); 
                   i_za_inc++)
                {
                  for (Index i_aa_inc = 0; i_aa_inc < pha_mat_data_raw.nrows(); 
                       i_aa_inc++)
                    {  
                      for (Index i = 0; i < pha_mat_data_raw.ncols(); i++)
                        {
                          pha_mat_data_int(i_za_sca, 
                                                i_aa_sca, i_za_inc, 
                                                i_aa_inc, i) =
                            interp(itw,
                                   pha_mat_data_raw(joker,i_za_sca, 
                                                i_aa_sca, i_za_inc, 
                                                i_aa_inc, i),
                                   freq_gp);
                        }
                    }
                }
            }
        }
      //       xml_write_to_file("pha_mat_data_unspr.xml", pha_mat_data);
      // xml_write_to_file("pha_mat_data_int.xml", pha_mat_data_int);
                          
       
       Matrix pha_mat_lab(stokes_dim, stokes_dim, 0.);
       
       
       // Do the transformation into the laboratory coordinate system.
       for (Index za_inc_idx = 0; za_inc_idx < scat_za_grid.nelem(); za_inc_idx ++)
         {
           for (Index aa_inc_idx = 0; aa_inc_idx < scat_aa_grid.nelem(); aa_inc_idx ++) 
             {
               pha_matTransform(pha_mat_lab,
                                pha_mat_data_int, 
                                za_datagrid, aa_datagrid,
                                part_type, scat_za_index, scat_aa_index, za_inc_idx, 
                                aa_inc_idx, scat_za_grid, scat_aa_grid,
                                scat_theta, scat_theta_gps, scat_theta_itws);
               pha_mat_spt(i_pt, za_inc_idx, aa_inc_idx, joker, joker) = pha_mat_lab;
             }
         }
    }
}
  

//! Calculates phase matrix for the single particle types.
/*! 
  In this function the phase matrix calculated
  in the laboratory frame. This function can be used in the agenda
  *pha_mat_spt_agenda* which is part of the calculation of the scattering
  integral (*scat_fieldCalc*).
 
  The function uses *pha_mat_sptDOITOpt*, which is the phase matrix interpolated
  on the right frequency for all scattering angles.
  
  The transformation from the database coordinate system to to 
  laboratory coordinate system is done in this function.

  NOTE: This is a special function for 1 particle type, including randomly 
  oriented particles.

  WS Output:
  \param pha_mat_spt Extinction matrix for a single particle type.
  WS Input:
  \param pha_mat_sptDOITOpt Phase matrix data.
  \param scat_theta Scattering angles.
  \param scat_za_grid Zenith angle grid for scattering calculation.
  \param scat_aa_grid Azimuth angle grid for scattering calculation.
  \param scat_za_index Index specifying the propagation direction.
  \param scat_aa_index  Index specifying the azimuth angle
  
  \author Claudia Emde
  \date   2003-08-25
*/

void pha_mat_sptFromDataDOITOpt( // Output:
                         Tensor5& pha_mat_spt,
                         // Input:
                         const Tensor6& pha_mat_sptDOITOpt,
                         const Tensor4& scat_theta,
                         const Index& za_grid_size,
                         const Vector& scat_aa_grid,
                         const Index& scat_za_index, // propagation directions
                         const Index& scat_aa_index
                         )
{

  Vector za_grid;
  nlinspace(za_grid, 0, 180, za_grid_size);  
  
  const Numeric za_sca = za_grid[scat_za_index];
  const Numeric aa_sca = scat_aa_grid[scat_aa_index];
  const Index N_pt = pha_mat_sptDOITOpt.nvitrines();

  Numeric za_inc;
  Numeric aa_inc;
  
  // Do the transformation into the laboratory coordinate system.
  for (Index i_pt = 0; i_pt < N_pt; i_pt ++)
    { 
      for (Index za_inc_idx = 0; za_inc_idx < za_grid_size; za_inc_idx ++)
        {
          for (Index aa_inc_idx = 0; aa_inc_idx < scat_aa_grid.nelem();
               aa_inc_idx ++) 
            {
              ConstVectorView pha_mat_int = pha_mat_sptDOITOpt
                (i_pt, scat_za_index, scat_aa_index, za_inc_idx, aa_inc_idx, joker);
              const Numeric theta_rad = scat_theta
                (scat_za_index, scat_aa_index, za_inc_idx, aa_inc_idx);
              
              za_inc = za_grid[za_inc_idx];
              aa_inc = scat_aa_grid[aa_inc_idx];
              
              pha_mat_labCalc(pha_mat_spt(i_pt, za_inc_idx, aa_inc_idx, joker, joker),
                              pha_mat_int, za_sca, aa_sca, za_inc, aa_inc, 
                              theta_rad);
            }
        }
    }

}


    
//! Calculates opticle properties for the single particle types.
/*! 
  In this function extinction matrix and absorption vector are calculated
  in the laboratory frame. These properties are required for the RT
  calculation, inside the the i_fieldUpdateXXX functions. 
  
  The interpolation of the data on the actual frequency is the first 
  step in this function. 
  
  Then the transformation from the database coordinate system to to 
  laboratory coordinate system is done.

  Output of the function are *ext_mat_spt*, and *abs_vec_spt* which hold 
  the optical properties for a specified propagation direction for 
  each particle type. 

  WS Output:
  \param ext_mat_spt Extinction matrix for a single particle type.
  \param abs_vec_spt Absorption vector for a single particle type.
  WS Input:
  \param scat_data_raw Raw data containing optical properties.
  \param scat_za_grid Zenith angle grid for scattering calculation.
  \param scat_aa_grid Azimuth angle grid for scattering calculation.
  \param scat_za_index Index specifying the propagation direction.
  \param scat_aa_index  Index specifying the azimuth angle
  \param f_index Index specifying the frequency.
  \param f_grid Frequency grid.

  \author Claudia Emde
  \date   2003-05-14
*/
void opt_prop_sptFromData( // Output and Input:
                         Tensor3& ext_mat_spt,
                         Matrix& abs_vec_spt,
                         // Input:
                         const ArrayOfSingleScatteringData& scat_data_raw,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid,
                         const Index& scat_za_index, // propagation directions
                         const Index& scat_aa_index,
                         const Index& f_index,
                         const Vector& f_grid
                         )
{
  const Index N_pt = scat_data_raw.nelem();
  const Index stokes_dim = ext_mat_spt.ncols();
  const Numeric za_sca = scat_za_grid[scat_za_index];
  const Numeric aa_sca = scat_aa_grid[scat_aa_index];

  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                         "must be 1,2,3 or 4");
  }
  
  assert( ext_mat_spt.npages() == N_pt );
  assert( abs_vec_spt.nrows() == N_pt );

  // Phase matrix in laboratory coordinate system. Dimensions:
  // [frequency, za_inc, aa_inc, stokes_dim, stokes_dim]
  Tensor3 ext_mat_data_int;
  Tensor3 abs_vec_data_int;
  
   // Initialisation
  ext_mat_spt = 0.;
  abs_vec_spt = 0.;


  // Loop over the included particle_types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      // First we have to transform the data from the coordinate system 
      // used in the database (depending on the kind of particle type 
      // specified by *ptype*) to the laboratory coordinate sytem. 
      
      // Frequency interpolation:
     
      // The data is interpolated on one frequency. 
      //
      // Resize the variables for the interpolated data:
      //
      ext_mat_data_int.resize(ext_mat_data_raw.npages(),
                              ext_mat_data_raw.nrows(), 
                              ext_mat_data_raw.ncols());
      //
      abs_vec_data_int.resize(abs_vec_data_raw.npages(),
                              abs_vec_data_raw.nrows(), 
                              abs_vec_data_raw.ncols());
      
      
      // Gridpositions:
      GridPos freq_gp;
      gridpos(freq_gp, f_datagrid, f_grid[f_index]); 

      // Interpolationweights:
      Vector itw(2);
      interpweights(itw, freq_gp);
     
      for (Index i_za_sca = 0; i_za_sca < ext_mat_data_raw.npages() ; i_za_sca++)
        {
          for(Index i_aa_sca = 0; i_aa_sca < ext_mat_data_raw.nrows(); i_aa_sca++)
            {
              //
              // Interpolation of extinction matrix:
              //
              for (Index i = 0; i < ext_mat_data_raw.ncols(); i++)
                {
                  ext_mat_data_int(i_za_sca, i_aa_sca, i) =
                    interp(itw, ext_mat_data_raw(joker, i_za_sca, i_aa_sca, i),
                           freq_gp);
                }
            }
        }

      for (Index i_za_sca = 0; i_za_sca < abs_vec_data_raw.npages() ; i_za_sca++)
        {
          for(Index i_aa_sca = 0; i_aa_sca < abs_vec_data_raw.nrows(); i_aa_sca++)
            {
              //
              // Interpolation of absorption vector:
              //
              for (Index i = 0; i < abs_vec_data_raw.ncols(); i++)
                {
                  abs_vec_data_int(i_za_sca, i_aa_sca, i) =
                    interp(itw, abs_vec_data_raw(joker, i_za_sca, i_aa_sca, i),
                           freq_gp);
                }
            }
        }
      

      //
      // Do the transformation into the laboratory coordinate system.
      //
      // Extinction matrix:
      //
     
  
      ext_matTransform(ext_mat_spt(i_pt, joker, joker),
                       ext_mat_data_int,
                       za_datagrid, aa_datagrid, part_type,
                       za_sca, aa_sca);
      // 
      // Absorption vector:
      //
      abs_vecTransform(abs_vec_spt(i_pt, joker),
                       abs_vec_data_int,
                       za_datagrid, aa_datagrid, part_type,
                       za_sca, aa_sca);                
    }


}
                          



//! Calculates absorption cross-section for a single particle

/*! 
  Calculates the absorption cross-section for a single particle 
  type from phase matrix and extinction  crosssection
  
  \param abs_vec_spt  Output : absorption vector for a single particle type
  \param ext_mat_spt  Input : extinction matrix for a single particle type
  \param pha_mat_spt  Input : phase matrix for  a single particle type
  \param scat_za_grid Input : zenith angle grid.
  \param scat_aa_grid Input : azimuth angle grid.
 
*/
void abs_vec_sptCalc(
                     Matrix& abs_vec_spt,
                     const Tensor3& ext_mat_spt,
                     const Tensor5& pha_mat_spt,
                     const Vector& scat_za_grid,
                     const Vector& scat_aa_grid)
                     
{
#ifndef NDEBUG
  Index nza = pha_mat_spt.nbooks(); 
  Index naa = pha_mat_spt.npages(); 
#endif
  Index npt = abs_vec_spt.nrows();
  Index stokes_dim = ext_mat_spt.nrows(); 

  assert (is_size(scat_za_grid, nza));
  assert (is_size(scat_aa_grid, naa));

  if ((stokes_dim > 4) || (stokes_dim <1)){
    throw runtime_error("The dimension of stokes vector "
                        "can be only 1,2,3, or 4");
  }
  
  if (abs_vec_spt.ncols() != stokes_dim ){
    //Dimension should agree to stokes_dim !!!
    throw runtime_error("The dimension of  abs_vec_spt should"
                        "agree to stokes_dim");
  }
  
  if ( ext_mat_spt.ncols() != stokes_dim ){
    throw runtime_error("The tensor ext_mat_spt should"
                        "  4 columns");
  }
  
  if (pha_mat_spt.nrows() != stokes_dim || 
      pha_mat_spt.ncols() != stokes_dim ){
    throw runtime_error("The tensor pha_mat_spt should"
                        "have 4 rows and 4 columns");
  }
  
  
  
  for (Index i = 0; i < npt; ++i)
    {
      ConstMatrixView ext = ext_mat_spt(i,
                                        Range(joker),
                                        Range(joker));
      
      ConstTensor4View pha=pha_mat_spt(i,
                                       Range(joker),
                                       Range(joker),
                                       Range(joker),
                                       Range(joker));
      
      // ConstVectorView za_grid = scat_za_grid[scat_za_index];
      //ConstVectorView aa_grid = scat_aa_grid[scat_aa_index];
      
      amp2abs(abs_vec_spt(i, Range(joker)),
              ext,
              pha,
              scat_za_grid,
              scat_aa_grid);
     
    }
}


//! Extinction Coefficient Matrix for the particle 
/*! 
  This function sums up the extinction matrices for all particle 
  types weighted with particle number density
  \param ext_mat Output and input : physical extinction coefficient 
  for the particles for given angles. 
  \param ext_mat_spt Input : extinction matrix for the single particle type
  \param pnd_field Input : particle number density givs the local 
  concentration for all particles.
  \param atmosphere_dim Input : he atmospheric dimensionality (now 1 or 3)
  \param scat_p_index Input : Pressure index for scattering calculations.
  \param scat_lat_index Input : Latitude index for scattering calculations.
  \param scat_lon_index Input : Longitude index for scattering calculations.
  
*/
void ext_matAddPart(
                      Tensor3& ext_mat,
                      const Tensor3& ext_mat_spt,
                      const Tensor4& pnd_field,
                      const Index& atmosphere_dim,
                      const Index& scat_p_index,
                      const Index& scat_lat_index,
                      const Index& scat_lon_index) 
                     
{
  Index N_pt = ext_mat_spt.npages();
  Index stokes_dim = ext_mat_spt.nrows();
  
  Matrix ext_mat_part(stokes_dim, stokes_dim, 0.0);

  
  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error(
                        "The dimension of stokes vector can be "
                        "only 1,2,3, or 4");
  }
  if ( ext_mat_spt.ncols() != stokes_dim){
    
    throw runtime_error(" The columns of ext_mat_spt should "
                        "agree to stokes_dim");
  }

  if (atmosphere_dim == 1)
    {
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt; l++)
        { 
          
          // now the last two loops over the stokes dimension.
          for (Index m = 0; m < stokes_dim; m++)
            {
              for (Index n = 0; n < stokes_dim; n++)
               //summation of the product of pnd_field and 
                //ext_mat_spt.
              ext_mat_part(m, n) += 
                (ext_mat_spt(l, m, n) * pnd_field(l, scat_p_index, 0, 0));
            }
        }

      //Add particle extinction matrix to *ext_mat*.
      ext_mat(0, Range(joker), Range(joker)) += ext_mat_part;
    }
 
  if (atmosphere_dim == 3)
    {
      
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt; l++)
        { 
          
          // now the last two loops over the stokes dimension.
          for (Index m = 0; m < stokes_dim; m++)
            {
              for (Index n = 0; n < stokes_dim; n++)
                 //summation of the product of pnd_field and 
                //ext_mat_spt.
                ext_mat_part(m, n) +=  (ext_mat_spt(l, m, n) * 
                                        pnd_field(l, scat_p_index, 
                                                  scat_lat_index, 
                                                  scat_lon_index));
              
            } 
        }

      //Add particle extinction matrix to *ext_mat*.
      ext_mat(0, Range(joker), Range(joker)) += ext_mat_part;

    }

} 

//! Absorption Vector for the particle 
/*! 
  This function sums up the absorption vectors for all particle 
  types weighted with particle number density
  \param abs_vec Output : physical absorption vector 
  for the particles for given angles. 
  \param abs_vec_spt Input : absorption for the single particle type
  \param pnd_field Input : particle number density givs the local 
  concentration for all particles.
  \param atmosphere_dim Input : the atmospheric dimensionality (now 1 or 3)
  \param scat_p_index Input : Pressure index for scattering calculations.
  \param scat_lat_index Input : Latitude index for scattering calculations.
  \param scat_lon_index Input : Longitude index for scattering calculations.
*/
void abs_vecAddPart(
                      Matrix& abs_vec,
                      const Matrix& abs_vec_spt,
                      const Tensor4& pnd_field,
                      const Index& atmosphere_dim,
                      const Index& scat_p_index,
                      const Index& scat_lat_index,
                      const Index& scat_lon_index) 
                    
{
  Index N_pt = abs_vec_spt.nrows();
  Index stokes_dim = abs_vec_spt.ncols();

  Vector abs_vec_part(stokes_dim, 0.0);

  if ((stokes_dim > 4) || (stokes_dim <1)){
    throw runtime_error("The dimension of stokes vector "
                        "can be only 1,2,3, or 4");
  } 
 
  if (atmosphere_dim == 1)
    {
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt ; ++l)
        {
          // now the loop over the stokes dimension.
          //(CE:) in the middle was l instead of m
          for (Index m = 0; m < stokes_dim; ++m)
             //summation of the product of pnd_field and 
            //abs_vec_spt.
            abs_vec_part[m] += 
              (abs_vec_spt(l, m) * pnd_field(l, scat_p_index, 0, 0));
          
        }
      //Add the particle absorption
      abs_vec(0, Range(joker)) += abs_vec_part;
    }
  
  if (atmosphere_dim == 3)
    {
      // this is a loop over the different particle types
      for (Index l = 0; l < N_pt ; ++l)
        {
          
          // now the loop over the stokes dimension.
          for (Index m = 0; m < stokes_dim; ++m)
             //summation of the product of pnd_field and 
            //abs_vec_spt.
            abs_vec_part[m] += (abs_vec_spt(l, m) *
                                pnd_field(l, scat_p_index,
                                          scat_lat_index, 
                                          scat_lon_index));
          
        }
      //Add the particle absorption
      abs_vec(0,Range(joker)) += abs_vec_part;
    }
} 

// //! Method for creating an absorption vector for gases for test purposes
// /*! 
  
// The aim of this method is just to give a reasonable value for absorption 
// coefficients that can be used to test the radiative transfer calculation.
// This method takes pressure grids and absorption values corresponding to them
// from ARTS.1.0 and interpolates it onto the pressure grid under consideration 
// for radiative transfer calculation.

// \param abs_vec_gas Output : The gaseous absorption vector
// \param p_grid Input : pressure grid on whiich the radiative transfer 
// calculations are done
// \param atmosphere_dim Input : atmospheric dimension(here  considered only 1D)
// \param stokes_dim Input : stokes dimension
// \param scat_p_index Input : The index corresponding to the pressure grid
// */
// void abs_vec_gasExample(Vector& abs_vec,
//                      const Vector& p_grid,
//                      const Index& atmosphere_dim,
//                      const Index& stokes_dim,
//                      const Index& scat_p_index) 
// {
//   abs_vec.resize( stokes_dim );
//   Vector abs_gas( p_grid.nelem() );
//   if (atmosphere_dim == 1){

//     //         MakeVector typical_abs(0,0,0,0,0,0,0,0,0,0);

//   //This is a typical absorption calculated from arts.1.0 (tropical)
  
//      // MakeVector typical_abs(0.0218978044006849,
// //                             0.0114605893154005,
// //                             0.00355854032448724,
// //                             0.00161793243125695,
// //                             0.000612486464168897,
// //                             0.000168382126027889,
// //                             3.35522552862195e-05,
// //                             8.52950996338938e-06,
// //                             4.32022347140681e-06,
// //                             3.12994113724593e-06);


//  //This is a typical absorption calculated from arts.1.0 (mls)
  
//     MakeVector typical_abs(0.0164163638663716,
//                         0.00768271579907592,
//                         0.00294668635075111,
//                         0.00125411825404778,
//                         0.000570445848162073,
//                         0.000236462958473072,
//                         4.40975932116215e-05,
//                         7.31218846316807e-06,
//                         3.643089167928e-06,
//                         3.12497184475723e-06);
    
// //This is a typical absorption calculated from arts.1.0 (mlw)
  
//     /*MakeVector typical_abs(0.00433659795310861,
//                         0.00267195350007522,
//                         0.00122676231151285,
//                         0.000524902357458338,
//                         0.000142742001685155,
//                         4.65799636740885e-05,
//                         1.07249251986862e-05,
//                         5.85896298848906e-06,
//                         4.94996692358466e-06,
//                         4.35302994835574e-06);*/


//     //The pressure grid for the above calculation (mlw)
//     MakeVector typical_abs_pgrid(100000,
//                               77426.3682681127,
//                               59948.4250318941,
//                               46415.8883361278,
//                               35938.1366380463,
//                               27825.5940220712,
//                               21544.3469003188,
//                               16681.0053720006,
//                               12915.4966501488,
//                               10000);
//     //p_grid is the new pressure grid
//     ArrayOfGridPos p_gp(p_grid.nelem());

//     gridpos(p_gp, typical_abs_pgrid, p_grid);

//     Matrix itw(p_grid.nelem(),2);

//     interpweights(itw, p_gp);
//     //interpolating absorption coefficients from typical_abs_pgrid to p_grid
//     interp(abs_gas, itw, typical_abs, p_gp);
      
//     for (Index i = 0; i<stokes_dim; ++i)
//       {
//      if (i == 0){
//        abs_vec[i] = abs_gas[scat_p_index];
//      }
//      else{
//        abs_vec[i] = 0.0;
//      }
//       }
//   }
  
// }

//! Initialize extinction matrix.
/*!
  This method is necessary, because all other extinction methods just
  add to the existing extinction matrix. 

  So, here we have to make it the right size and fill it with 0.
  
  Note, that the matrix is not really a matrix, because it has a
  leading frequency dimension.

  \param ext_mat Extinction matrix.
  \param f_grid Frequency grid.
  \param stokes_dim Stokes dimension.
  \param f_index Frequency index.

  \author Stefan Buehler

  \date   2002-12-12
*/
void ext_matInit( Tensor3&      ext_mat,
                  const Vector& f_grid,
                  const Index&  stokes_dim,
                  const Index&   f_index)
{
  Index freq_dim;

  if( f_index < 0 )
    freq_dim = f_grid.nelem();
  else
    freq_dim = 1;
 
  ext_mat.resize( freq_dim,
                  stokes_dim,
                  stokes_dim );
  ext_mat = 0;                  // Initialize to zero!

  out2 << "Set dimensions of ext_mat as ["
       << freq_dim << ","
       << stokes_dim << ","
       << stokes_dim << "] and initialized to 0.\n";
}

//! Add gas absorption to all diagonal elements of extinction matrix.
/*! 
  The task of this method is to sum up the gas absorption of the
  different gas species and add the result to the extinction matrix. 

  \param ext_mat Input and Output: Extinction matrix.
                 Dimension: [#frequencies, stokes_dim, stokes_dim]
  
  \param abs_scalar_gas Scalar gas absorption coefficients.
                        Dimension: [#frequencies, #gas_species]

  \author Stefan Buehler
          (rewritten version of function by Sreerekha Ravi)

  \date   2002-12-12
*/
void ext_matAddGas( Tensor3&      ext_mat,
                    const Matrix& abs_scalar_gas )
{
  // Number of Stokes parameters:
  const Index stokes_dim = ext_mat.ncols();

  // The second dimension of ext_mat must also match the number of
  // Stokes parameters:
  if ( stokes_dim != ext_mat.nrows() )
    throw runtime_error("Row dimension of ext_mat inconsistent with "
                        "column dimension."); 

  // Number of frequencies:
  const Index f_dim = ext_mat.npages();

  // This must be consistent with the first dimension of
  // abs_scalar_gas. Check this:
  if ( f_dim != abs_scalar_gas.nrows() )
    throw runtime_error("Frequency dimension of ext_mat and abs_scalar_gas\n"
                        "are inconsistent.");

  // Sum up absorption over all species.
  // This gives us an absorption vector for all frequencies. Of course
  // this includes the special case that there is only one frequency.
  Vector abs_total(f_dim);
  for ( Index i=0; i<f_dim; ++i )
    abs_total[i] = abs_scalar_gas(i,joker).sum();

  for ( Index i=0; i<stokes_dim; ++i )
    {
      // Add the absorption value to all the diagonal elements:
      ext_mat(joker,i,i) += abs_total[joker];
      
      // We don't have to do anything about the off-diagonal
      // elements! 
    }
}

//! Initialize absorption vector.
/*!
  This method is necessary, because all other absorption methods just
  add to the existing absorption vector.

  So, here we have to make it the right size and fill it with 0.

  Note, that the vector is not really a vector, because it has a
  leading frequency dimension.
  
  \param abs_vec Extinction matrix.
  \param f_grid Frequency grid.
  \param stokes_dim Stokes dimension.
  \param f_index Frequency index.

  \author Stefan Buehler

  \date   2002-12-12
*/
void abs_vecInit( Matrix&       abs_vec,
                  const Vector& f_grid,
                  const Index&  stokes_dim,
                  const Index&  f_index)
{
  Index freq_dim;

  if( f_index < 0 )
    freq_dim = f_grid.nelem();
  else
    freq_dim = 1;
 
  abs_vec.resize( freq_dim,
                  stokes_dim );
  abs_vec = 0;                  // Initialize to zero!

  out2 << "Set dimensions of abs_vec as ["
       << freq_dim << ","
       << stokes_dim << "] and initialized to 0.\n";
}

//! Add gas absorption to first element of absorption vector.
/*! 
  The task of this method is to sum up the gas absorption of the
  different gas species and add the result to the first element of the
  absorption vector.

  \param abs_vec Input and Output: Absorption vector
                 Dimension: [#frequencies, stokes_dim]
  
  \param abs_scalar_gas Scalar gas absorption coefficients.
                        Dimension: [#frequencies, #gas_species]

  \author Stefan Buehler
          (rewritten version of function by Sreerekha Ravi)

  \date   2002-12-12
*/
void abs_vecAddGas( Matrix&       abs_vec,
                    const Matrix& abs_scalar_gas )
{
  // Number of frequencies:
  const Index f_dim = abs_vec.nrows();

  // This must be consistent with the first dimension of
  // abs_scalar_gas. Check this:
  if ( f_dim != abs_scalar_gas.nrows() )
    throw runtime_error("Frequency dimension of abs_vec and abs_scalar_gas\n"
                        "are inconsistent.");

  // Loop all frequencies. Of course this includes the special case
  // that there is only one frequency.
  for ( Index i=0; i<f_dim; ++i )
    {
      // Sum up the columns of abs_scalar_gas and add to the first
      // element of abs_vec.
      abs_vec(i,0) += abs_scalar_gas(i,joker).sum();
    }

  // We don't have to do anything about higher elements of abs_vec,
  // since scalar gas absorption only influences the first element.
}

//! Add Zeeman extinction  to the elements of extinction matrix.
/*! 
  
\param ext_mat Input and Output: Extinction matrix.
Dimension: [#frequencies, stokes_dim, stokes_dim]

\param ext_mat_zee Input : Zeeman extinction coefficient matrix.
Dimension: [#frequencies, stokes_dim, stokes_dim]

\author Sreerekha Ravu
          (based on ext_matAddGas by Stefan Buehler)

  \date   2003-12-01
*/
void ext_matAddGasZeeman( Tensor3&      ext_mat,
			  const Tensor3&  ext_mat_zee )
{
  // Number of Stokes parameters:
  const Index stokes_dim = ext_mat.ncols();

  // The second dimension of ext_mat must also match the number of
  // Stokes parameters:
  if ( stokes_dim != ext_mat.nrows() )
    throw runtime_error("Row dimension of ext_mat inconsistent with "
                        "column dimension."); 

  for ( Index i=0; i<stokes_dim; ++i )
    {
      for ( Index j=0; j<stokes_dim; ++j )
	{
	  // Add the zeeman extinction to extinction matrix.
	  ext_mat(joker,i,j) += ext_mat_zee(joker, i, j);
	}
      
    }
}

//! Add zeeman absorption to the elements of absorption vector.
/*! 
  The task of this method is to sum up the gas absorption of the
  different gas species and add the result to the first element of the
  absorption vector.

  \param abs_vec Input and Output: Absorption vector
                 Dimension: [#frequencies, stokes_dim]
  
  \param abs_vec_zee Input : Zeeman absorption vector.
                        Dimension: [#frequencies, #gas_species]

  \author Sreerekha Ravi
          (based on abs_vecAddGas written by Stefan Buehler)

  \date   2003-12-01
*/
void abs_vecAddGasZeeman( Matrix&      abs_vec,
			  const Matrix& abs_vec_zee )
{
  // Number of Stokes parameters:
  const Index stokes_dim = abs_vec_zee.ncols();
  // that there is only one frequency.
  for ( Index j=0; j<stokes_dim; ++j )
    {
      abs_vec(joker,j) += abs_vec_zee(joker,j);
    }
}


//! Phase Matrix for the particle 
/*! 
  
  This function sums up the phase matrices for all particle 
  types weighted with particle number density
  \param pha_mat Output : physical phase matrix 
  for the particles for given angles. 
  \param pha_mat_spt Input : phase matrix for the single particle type
  \param pnd_field Input : particle number density givs the local 
  concentration for all particles.
  \param atmosphere_dim Input : he atmospheric dimensionality (now 1 or 3)
  \param scat_p_index Input : Pressure index for scattering calculations.
  \param scat_lat_index Input : Latitude index for scattering calculations.
  \param scat_lon_index Input : Longitude index for scattering calculations.
*/

void pha_matCalc(
                 Tensor4& pha_mat,
                 const Tensor5& pha_mat_spt,
                 const Tensor4& pnd_field,
                 const Index& atmosphere_dim,
                 const Index& scat_p_index,
                 const Index& scat_lat_index,
                 const Index& scat_lon_index) 
                      
{

  Index N_pt = pha_mat_spt.nshelves();
  Index Nza = pha_mat_spt.nbooks();
  Index Naa = pha_mat_spt.npages();
  Index stokes_dim = pha_mat_spt.nrows();
 
  pha_mat.resize(Nza, Naa, stokes_dim, stokes_dim);

  // Initialisation
  for (Index za_index = 0; za_index < Nza; ++ za_index)
    {
      for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
        {
          for (Index stokes_index_1 = 0; stokes_index_1 < stokes_dim;
               ++ stokes_index_1)
            {
              for (Index stokes_index_2 = 0; stokes_index_2 < stokes_dim; 
                   ++ stokes_index_2)
                pha_mat(za_index, aa_index, stokes_index_1, stokes_index_2)

                  = 0.0;
            }
        }
    }
  if (atmosphere_dim == 1)
    {
      // this is a loop over the different particle types
      for (Index pt_index = 0; pt_index < N_pt; ++ pt_index)
        {
                  
          // these are loops over zenith angle and azimuth angle
          for (Index za_index = 0; za_index < Nza; ++ za_index)
            {
              for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
                {
                  
                  // now the last two loops over the stokes dimension.
                  for (Index stokes_index_1 = 0; stokes_index_1 < stokes_dim; 
                       ++  stokes_index_1)
                    {
                      for (Index stokes_index_2 = 0; stokes_index_2 < stokes_dim;
                           ++ stokes_index_2)
                         //summation of the product of pnd_field and 
                          //pha_mat_spt.
                        pha_mat(za_index, aa_index,  
                                     stokes_index_1, stokes_index_2) += 
                          
                          (pha_mat_spt(pt_index, za_index, aa_index,  
                                       stokes_index_1, stokes_index_2) * 
                           pnd_field(pt_index,scat_p_index, 0, 0));
                    }
                }
            }
        }
    }
          
  if (atmosphere_dim == 3)
    {
      // this is a loop over the different particle types
      for (Index pt_index = 0; pt_index < N_pt; ++ pt_index)
        {
          
          // these are loops over zenith angle and azimuth angle
          for (Index za_index = 0; za_index < Nza; ++ za_index)
            {
              for (Index aa_index = 0; aa_index < Naa; ++ aa_index)
                {
                  
                  // now the last two loops over the stokes dimension.
                  for (Index i = 0;  i < stokes_dim; ++  i)
                    {
                      for (Index j = 0; j < stokes_dim; ++ j)
                        {
                          //summation of the product of pnd_field and 
                          //pha_mat_spt.
                          pha_mat(za_index, aa_index, i,j ) += 
                            (pha_mat_spt(pt_index, za_index, aa_index, i, j) * 
                             pnd_field(pt_index, scat_p_index,  
                                       scat_lat_index, scat_lon_index));
                          
                          
                        } 
                    }	
                }
            }		
        }	
    }		
}





// //! Extinction coefficient matrix  for convergence test.
// /*! 
//   This function sums up the convergence extinction matrices for all particle 
//   types weighted with particle number density. This is exactly the same 
//   function as *ext_mat_partCalc*, only with different input, *ext_mat_conv_spt*
//   instead of ext_mat_spt.

//   \param ext_mat_part Output : physical extinction coefficient 
//   for the particles for given angles. 
//   \param ext_mat_conv_spt Input : extinction matrix for the single 
//   particle type
//   \param pnd_field Input : particle number density givs the local 
//   concentration for all particles.
//   \param atmosphere_dim Input : he atmospheric dimensionality (now 1 or 3)
//   \param scat_p_index Input : Pressure index for scattering calculations.
//   \param scat_lat_index Input : Latitude index for scattering calculations.
//   \param scat_lon_index Input : Longitude index for scattering calculations.
// */
// void ext_mat_AddScat(
//                    Matrix& ext_mat,
//                    const Tensor3& ext_mat_spt,
//                    const Tensor4& pnd_field,
//                    const Index& atmosphere_dim,
//                    const Index& scat_p_index,
//                    const Index& scat_lat_index,
//                    const Index& scat_lon_index) 
                     
// {
//   Index N_pt = ext_mat_spt.npages();
//   Index stokes_dim = ext_mat_spt.nrows();
  
//   Matrix ext_mat_part(stokes_dim, stokes_dim, 0.);
  
 
//   if (atmosphere_dim == 1)
//     {
//       // this is a loop over the different particle types
//       for (Index l = 0; l < N_pt; l++)
//         { 
       
//        // now the last two loops over the stokes dimension.
//        for (Index m = 0; m < stokes_dim; m++)
//          {
//            for (Index n = 0; n < stokes_dim; n++)
//              //summation of the product of pnd_field and 
//              //ext_mat_conv_spt.
//              ext_mat_part(m, n) += 
//              (ext_mat_spt(l, m, n) * pnd_field(l, scat_p_index, 0, 0));
//          }
//      }
//       //Add particle extinction matrix to *ext_mat*.
//       ext_mat(0, Range(joker), Range(joker)) += ext_mat_part;
//     }
//   if (atmosphere_dim == 3)
//     {
      
//       // this is a loop over the different particle types
//       for (Index l = 0; l < N_pt; l++)
//      { 
          
//        // now the last two loops over the stokes dimension.
//        for (Index m = 0; m < stokes_dim; m++)
//          {
//            for (Index n = 0; n < stokes_dim; n++)
//               //summation of the product of pnd_field and 
//              //ext_mat_conv_spt.
//              ext_mat_part(m, n) +=  (ext_mat_spt(l, m, n) * 
//                                      pnd_field(l, scat_p_index, 
//                                                scat_lat_index, 
//                                                scat_lon_index));
              
//          } 
//      }
//     }
// } 


//! Method for getting amp_mat from amp_mat_raw.
/*! 
  This method converts the raw amplitude matrix data namely amp_mat_raw 
  to workspace variable amp_mat which can be directly used for calculating  
  scattering properties of particles.  The data type of amp_mat_raw is 
  an ArrayOfArrayOfTensor6 which contains one gridded field for each 
  particle type the user wants.  One data file
  contains amp_mat_raw for one particle type.  One data file has
  amp_mat_raw for one particle for a set of frequencies, incoming and 
  outgoing zenith angles and azimuth angles. The frequencies, angles, 
  and the amplitude matrix are each a Tensor 6.  The size of amp_mat_raw 
  is amp_mat_raw[Npt][7].  where Npt is the number of particle types.
  amp_mat_raw[Npt][0] gives the frequency tensor [N_f, 1, 1, 1, 1, 1] 
  where N_f gives the number of frequencies considered in that particular
  database file. Similarly,
  amp_mat_raw[Npt][1] gives the outgoing zenith angle tensor [1,Nza,1,1,1,1],
  amp_mat_raw[Npt][2] gives the outgoing azimuth angle tensor [1,1,Naa,1,1,1], 
  amp_mat_raw[Npt][3] gives the incoming zentih angle tensor [1,1,1,Nza,1,1], 
  amp_mat_raw[Npt][4] gives the incoming azimuth angle tensor [1,1,1,1,Naa,1],
  amp_mat_raw[Npt][5] is a dummy tensor6 and 
  amp_mat_raw[Npt][6] gives amplitude matrix which is also a tensor 6 of
  size [N_f, N_za, N_aa, N_za, N_aa, 8]. Here, Nza is the number 
  of zenith angles, Naa is the number of azimuth angles and 8 denotes the 
  amplitude matrix elements.  
  
  In this method, we have to interpolate the raw data calculated on a 
  specific angular and frequency grid onto a grid which is set by the user. 
  Since we decide that frequency should be the outermost loop for our 
  radiative transfer calculation the frequency grid contains just one 
  value specified by the index f_index.  The angles for which the 
  calculations are to be done are specified by scat_za_grid and scat_aa_grid.
  
  The output of this method is amp_mat has to be a 
  Tensor6 with the first dimension being that of the particle type, then the
  angles and finally the amplitude matrix element 8. The size of amp_mat is 
  (Npt, Nza, Naa, Nza, Naa, 8).  Note that the dimension frequency is taken
  out.
  
  \param amp_mat Output : Amplitude matrix which is interpolated on the 
  required grids set by the user.
  \param amp_mat_raw Input : The original data as read from the database.
  \param f_index Input : The frequency index.
  \param f_grid Input : The frequency grid as required for RT calculation
  \param scat_za_grid Input : The zenith angle grid for which scattering 
  properties are to be calculated 
  \param scat_aa_grid Input : The azimuth angle grid for which scattering 
  properties are to be calculated 
*/
void amp_matCalc(Tensor6& amp_mat,
                 const ArrayOfArrayOfTensor6& amp_mat_raw,
                 const Index& f_index,
                 const Vector& f_grid,
                 const Vector& scat_za_grid,
                 const Vector& scat_aa_grid)
{
  Index N_pt = amp_mat_raw.nelem();
  
  Index N_za = scat_za_grid.nelem();
  Index N_aa = scat_aa_grid.nelem();
  Index N_0 = amp_mat_raw [ 0 ] [ 6 ].ncols();
  Index nza1_from = amp_mat_raw [ 0 ]  [ 1 ].nshelves();
  
 
  if (N_0 != 8)
    throw runtime_error(
                        "Amplitude matrix must have 8 columns because there"
                        "are 4 real and 4 complex components stored"
                        "separately");

  amp_mat.resize(N_pt, N_za, N_aa, N_za, N_aa, N_0);  
  /*amp_mat_raw is an ArrayOfArrayOfTensor6 from which we can get the 
    amplitude matrix which is a tensor 6 with size(Nf, Nza, Naa, Nza,Naa,8). 
    This is what is directly read from the database.  As mentioned in 
    the beginning in this method we have to interpolate amplitude matrix 
    from the grid on which it is calculated in the database to the grid set
    by the user. The interpolation has to be done for the frequency grid,
    zenith angle grid and azimuth angle grid. Since this interpolation is
    from a gridded field to a new field, it is called a green interpolation.
    For more insight into the interpolation schemes refer to Chapter 8-
    Interpolation  of AUG.
    The new grids are f_grid(only one element in this vector), scat_za_grid, 
    scat_aa_grid and the old grids are deroved from the ArrayOfArrayOfTensor6 
    amp_amt_raw.*/

  //Loop over the particle types. We can get information about the number of
  //particle types from the input ArrayOfArrayOfTensor6 amp_mat_raw
 
 
  if (is_size(amp_mat_raw[0][2],
	      0,0,0,0,0,0) &&
      is_size(amp_mat_raw[0][3],
	      0,0,0,0,0,0) &&
      is_size(amp_mat_raw[0][4],
	      0,0,0,0,0,0)) 
    {
      
      for (Index ipt = 0; ipt < N_pt; ++ ipt )
	{
	  
	  Index N_i = amp_mat_raw [ ipt ] [ 6 ].ncols();
	  //Tensor4 za_big(N_za, N_aa, N_za, N_aa);
	  //Vector za_big(N_aa * N_za *N_aa *N_za);
	  Numeric za_big = 0.0;
	  cout <<  " nza1_from" <<  nza1_from<<endl;
	  //calling the interpolation routines.  
	  for (Index za_out = 0; za_out < N_za  ; ++ za_out) 
	    {
	      for (Index aa_out = 0; aa_out < N_aa  ; ++ aa_out)
		{
		  for (Index za_in = 0; za_in < N_za ; ++ za_in)
		    {
		      for (Index aa_in = 0; aa_in < N_aa ; ++ aa_in)
			{
			  za_big = 
			    acos(
				 (cos (scat_za_grid[za_out] * PI/180.) * 
				  cos (scat_za_grid[za_in] * PI/180.)) +
				 (sin (scat_za_grid[za_out] * PI/180.) *
				  sin (scat_za_grid[za_in] * PI/180.) *
				  cos((scat_aa_grid[aa_in] * PI/180.) - 
				      (scat_aa_grid[aa_out] * PI/180.)))) *
			    180./PI; 
			  
			 		
			  // cout << "  za _ big" << "  " << za_big << endl;
			  
			  // for  frequency : 
			  ArrayOfGridPos freq_gp(1); 
			  
			  gridpos (freq_gp,
				   amp_mat_raw [ 0 ]  [ 0 ] 
				   ( Range(joker), 0, 0, 0, 0, 0),
				   f_grid[Range(f_index, 1)]);
			  
			  //like for frquency we can get the 
			  //gridpostions for the angular grids.
 			  ArrayOfGridPos za_big_gp(1);
 			  gridpos (za_big_gp, 
 				   amp_mat_raw [ ipt ] [  1 ]
				   ( 0, Range(joker), 0, 0, 0, 0),
 				   za_big);
			  // 		
 			  Matrix itw(1,2);
			  
			  interpweights ( itw, za_big_gp);

 			  Tensor6 target_amp(N_pt,
					     N_za * N_za * N_aa * N_aa,
					     1,
					     1,
					     1,
					     8);

 			  for (Index i = 0 ; i < N_i ; ++ i)
			    {
			      
			      VectorView target_amp_view = 
				target_amp(Range(ipt, 1),
					   1,
					   0,
					   0,
					   0,
					   i);
			      
			      ConstVectorView source_amp = 
				amp_mat_raw [ ipt ][6](f_index,
						       Range(joker),
						       0,
						       0,
						       0,
						       i); 
			      interp (target_amp_view, 
				      itw, 
				      source_amp, 
				      za_big_gp );
			      amp_mat(Range(joker),
				      za_out,
				      aa_out,
				      za_in,
				      aa_in,
				      i) = target_amp_view;
			      //cout << target_amp_view<< endl;
			      //}// closes 2nd j loop
			    }// closes i loop
			}// closes aa_in loop
		    }// closes za_in loop
		}// closes aa_out loop
	    }// closes za_out loop 		
	  
	}// closes particle type loop
    }// closes if condition for macroscopically symmetric and isotropic particles
  else
    {
      for (Index ipt = 0; ipt < N_pt; ++ ipt )
	{
	  Index N_i = amp_mat_raw [ ipt ] [ 6 ].ncols();
	  
	  //Define the grid position arrays. 
	  
	  // for  frequency : 
	  ArrayOfGridPos f_gp(1); 
	  
	  // for outgoing zenith angle grids : 
	  ArrayOfGridPos za_gp(scat_za_grid.nelem());
	  
	  // for outgoing azimuth angle grids : 
	  ArrayOfGridPos aa_gp(scat_aa_grid.nelem());
	  
	  // for incoming zenith angle grids : 
	  ArrayOfGridPos za_in_gp(scat_za_grid.nelem());
	  
	  // for incoming azimuth angle grids : rows
	  ArrayOfGridPos aa_in_gp(scat_aa_grid.nelem());
	  
	  // Set up Grid position arrays by calling the function gridpos.
	  
	  /*for frequency :
	    
	  f_gp is the ArrayOfGridpos.
	  
	  original frequency grid as in the data base can be got from the 
	  ArrayOfTensor6 amp_mat_raw.  
	  amp_mat_raw[ipt][0] ( Range(joker), 0, 0, 0, 0, 0).  
	  
	  f_grid(Range(f_index),1) is the new grid which is also a vector
	  formally but with one element given by the f_index*/
	  //cout<<f_index<<"\n";
	  gridpos (f_gp,
		   amp_mat_raw [ ipt ]  [ 0 ] ( Range(joker), 0, 0, 0, 0, 0),
		   f_grid[Range(f_index, 1)]);
	   
	  //like for frquency we can get the gridpostions for the angular grids.
	  
	  gridpos (za_gp, 
		   amp_mat_raw [ ipt ] [  1 ] ( 0, Range(joker), 0, 0, 0, 0),
		   scat_za_grid);
	  
	  gridpos (aa_gp,
		   amp_mat_raw [ ipt ] [ 2 ] ( 0, 0, Range(joker), 0, 0, 0),
		   scat_aa_grid);
	  
	  gridpos (za_in_gp, 
		   amp_mat_raw [ ipt ] [ 3 ] ( 0, 0, 0, Range(joker), 0, 0),
		   scat_za_grid);
	  
	  gridpos (aa_in_gp,
		   amp_mat_raw [ ipt ] [ 4] ( 0, 0, 0, 0, Range(joker), 0),
		   scat_aa_grid);
	  
	  /*The interpolation weights. Since here we interpolate 
	    simultaneously in 5 dimensions we require exactly 32 
	    interpolation weights. (There are 2^n interpolation 
	    weights for an n-dimensional interpolation. ).  Since
	    in ARTS we want to save a lot by re-using the weights
	    the interpolation weights are stored in a tensor whcih
	    has one more dimension than the output field. The last 
	    dimension is for weight, this explains the last dimension
	    of interpolation weight to be 32 in this case.
	  */
	  
	  Tensor6 itw(1, scat_za_grid.nelem(), scat_aa_grid.nelem(),
		      scat_za_grid.nelem(), scat_aa_grid.nelem(), 32);
	  
	  //function for computing interpolation weight tensor. For this
	  //step also we need only gridpositions, not the field yet.
	  interpweights ( itw, f_gp, za_gp, aa_gp, za_in_gp, aa_in_gp );
	  
	  // This is a green interpolation for all columns of source_amp.
	  // Loop over the column
	  for (Index i = 0 ; i < N_i ; ++ i)
	    {
	      /*Select the current column of target.  The first element 
		of target is the running variable ipt, whcih gives the
		particle type.*/
	      Tensor5View target_amp = amp_mat(Range(ipt, 1),
					       Range(joker),
					       Range(joker),
					       Range(joker),
					       Range(joker),
					       i);
	      
	      //Select the current column of source.  Here the first
	      //element is frequency
	      ConstTensor5View source_amp = 
		amp_mat_raw [ ipt ][6 + 7* ipt](Range(joker),
						Range(joker),
						Range(joker),
						Range(joker),
						Range(joker),
						i);
	      
	      //Interpolation is done :
	      interp (target_amp, 
		      itw, 
		      source_amp, 
		      f_gp, 
		      za_gp,
		      aa_gp,
		      za_in_gp,
		      aa_in_gp);
	      
	    }//close column index
	}//close particle index loop
      //cout<< amp_mat<<"\n";
    }
}


//! Check the scattering data in the database.
/*! 
  This function can be used to check datafiles containing data for 
  randomly oriented scattering media.
  It is checked whether the data is consistent. The integral over the phase 
  matrix should result the scattering cross section \<C_sca\>.
  
  The check is if:
  \<C_ext\> - \<C_sca\> = \<C_abs\>
  
  The result is printed on the screen.
  
  \param scat_data_raw (Input) Single scattering data.
  
  \author Claudia Emde
  \date   2003-05-20
*/
void scat_data_rawCheck(//Input:
                         const ArrayOfSingleScatteringData& scat_data_raw
                         )
{

  xml_write_to_file("SingleScatteringData", scat_data_raw);
  
  const Index N_pt = scat_data_raw.nelem();
  
  // Loop over the included particle_types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      Numeric Csca = AngIntegrate_trapezoid
        (pha_mat_data_raw(0, joker, 0, 0, 0, 0), za_datagrid);

      Numeric Cext = ext_mat_data_raw(0,0,0,0);

      Numeric Cabs = Cext - Csca;

      Numeric Cabs_data = abs_vec_data_raw(0,0,0,0);

      Numeric Csca_data = Cext - Cabs_data;

      out1 << " Coefficients in database: \n"
           << "Cext: " << Cext << " Cabs: " << Cabs_data << " Csca: " << Csca_data  
           << " \n Calculated absorption cooefficient: \n"
           << "Cabs calculated: " << Cabs   
           << " Csca: " << Csca << "\n";
      
    }
}

//! Prepare single scattering data for a DOIT scattering calculation.
/*! 
  This function has to be used for scattering calculations using the 
  DOIT method. It prepares the data in such a way that the code can be 
  speed-optimized.
  
  For different hydrometeor species different preparations are required. 
  So far, we only the case of randomly oriented particles is implemented. 
  All files of hydrometeor must have the same angular grids. Otherwise
  the output variables of this method would become to large.

  For randomly oriented hydrometeor species the scattering angles in the 
  scattering frame are precalculated for all possible incoming and 
  scattered directions and stored the the WSV *scat_theta*. In the program 
  all phase matrix emements have to be interpolated on the scattering angle.
  This has to be done for each iteration and for all frequencies.
  
  Output:
  \param scat_theta All scattering angles.
  \param scat_theta_gps Grid positions of scattering angles.
  \param scat_theta_itws Interpolation weights.
  Input:
  \param scat_za_grid Zenith angle grid for scattering calculations.
  \param scat_aa_grid Azimuth angle grid for scattering calculations.
  \param scat_data_raw Array of single scattering data for all hydrometeor
                      species

  \author Claudia Emde
  \date   2003-08-20
*/
void ScatteringDataPrepareDOIT( //Output:
                               Tensor4& scat_theta,
                               ArrayOfArrayOfArrayOfArrayOfGridPos& 
                               scat_theta_gps,
                               Tensor5& scat_theta_itws,
                               //Input:
                               const Vector& scat_za_grid,
                               const Vector& scat_aa_grid,
                               const ArrayOfSingleScatteringData& scat_data_raw
                    )
{
  // At the moment this function can only be used for randomly orientes 
  // particles.
  // FIXME: Include other cases!!!
 
  // Check, whether the angular grids of all data types are the same:
  // for ( Index i_pt = 1; i_pt < scat_data_raw.nelem(); i_pt ++)
//     {
//       if( scat_data_raw[i_pt].za_grid != scat_data_raw[i_pt-1].za_grid)
//         throw runtime_error("The angular grids of all hydrometeor species\n"
//                             "must be identical to use this method. \n"
//                             "If not you can use *ScatteringDataPrepareOff*\n"
//                             "which works for all file. \n"
//                             );
//     }
     
  // Initialize variables:
  scat_theta.resize(scat_za_grid.nelem(),scat_aa_grid.nelem(), 
                    scat_za_grid.nelem(),scat_aa_grid.nelem());
  
  scat_theta_gps.resize(scat_za_grid.nelem());

  for( Index i = 0; i<scat_za_grid.nelem(); i++)
    scat_theta_gps[i].resize(scat_aa_grid.nelem());
  
  for( Index i = 0; i<scat_za_grid.nelem(); i++)
    {
      for( Index j = 0; j<scat_aa_grid.nelem(); j++)
        {
          scat_theta_gps[i][j].resize(scat_za_grid.nelem());
        }
    }
     
  for( Index i = 0; i<scat_za_grid.nelem(); i++)
    {
      for( Index j = 0; j<scat_aa_grid.nelem(); j++)
        {
          for( Index k = 0; k<scat_za_grid.nelem(); k++)
            {
              scat_theta_gps[i][j][k].resize(scat_aa_grid.nelem());
            }
        }
    }

  scat_theta_itws.resize(scat_za_grid.nelem(),scat_aa_grid.nelem(), 
                         scat_za_grid.nelem(),scat_aa_grid.nelem(),2);
  

  // Calculate all scattering angles for all combinations of incoming 
  // and scattered directions.
  
  // To calculate the gridpositions it is more efficient 
  // to calculate them for all scattering angles at the 
  // same time. For this, we need the scattering angles in 
  // a sorted vector.
  
  // First partcle type is taken.
  Index i_pt = 0;

  // First we create a vector including all scattering angles:
  Index N_stv = scat_za_grid.nelem() * scat_za_grid.nelem() *
    scat_aa_grid.nelem() * scat_aa_grid.nelem();
  Vector scat_theta_vector(N_stv);
  
  // Counting index for elements of vector containing scattering angles.
  Index i = 0;
  
  Numeric za_sca_rad, za_inc_rad, aa_sca_rad, aa_inc_rad;
  Numeric theta_rad;
  
  cout << "Calculate scattering angles:" << endl;
  
   for (Index za_sca = 0; za_sca < scat_za_grid.nelem(); za_sca ++)
    {
      for (Index aa_sca = 0; aa_sca < scat_aa_grid.nelem(); aa_sca ++)
        {
          for (Index za_inc = 0; za_inc < scat_za_grid.nelem(); za_inc ++)
            {
              for (Index aa_inc = 0; aa_inc < scat_aa_grid.nelem(); aa_inc ++)
                {
                  za_sca_rad = scat_za_grid[za_sca] * DEG2RAD;
                  za_inc_rad = scat_za_grid[za_inc] * DEG2RAD;
                  aa_sca_rad = scat_aa_grid[aa_sca] * DEG2RAD;
                  aa_inc_rad = scat_aa_grid[aa_inc] * DEG2RAD;
                  
                  theta_rad = acos(cos(za_sca_rad) * cos(za_inc_rad) + 
                                   sin(za_sca_rad) * sin(za_inc_rad) * 
                                   cos(aa_sca_rad - aa_inc_rad));
                  
                  scat_theta(za_sca, aa_sca, za_inc, aa_inc)
                    = theta_rad;
                  
                  // cout << "scat_theta: " << scat_theta(za_sca, aa_sca, za_inc, aa_inc) << endl;
                  
                  scat_theta_vector[i] = RAD2DEG * theta_rad;
                  
                  // Increase i 
                  i++;
                }
            }
        }
    }
    
   // Sort scat_theta_vector.

   cout << "Sort the scattering angles" << endl;

   // i_sorted is a vector containing the sorting indices.
   ArrayOfIndex i_sorted;
   get_sorted_indexes(i_sorted, scat_theta_vector);

   Vector sorted_scat_theta_vector(N_stv);

   for (Index j = 0; j < N_stv; j++)
     sorted_scat_theta_vector[i_sorted[j]] = scat_theta_vector[j];
            
            

   cout<< "Initializae gps, itw" << endl;
     
   ArrayOfGridPos gps_sorted_array(N_stv);
   Matrix itws_sorted_array(N_stv, 2);
   
   cout << "calculate gridpos" << endl;
   
   gridpos(gps_sorted_array, za_datagrid, sorted_scat_theta_vector);

   cout << "calculate itw" << endl;
   interpweights(itws_sorted_array, gps_sorted_array);
   
   cout << "initialize" << endl;
   // Sort gridpostions and interpolation weights back into the original order:
   ArrayOfGridPos gps_array(N_stv);
   Matrix itws_array(N_stv, 2);
   
   
   
   cout << "back- sorting" << endl;

   for (Index j = 0; j < N_stv; j++)
     {
       gps_array[i_sorted[j]] = gps_sorted_array[j];
       itws_array(i_sorted[j], joker) = itws_sorted_array(j, joker);
       //       cout << " sorted array " << gps_sorted_array[j] << endl;
     }
 
   
   // Put the grid positions and interpolation weights into the WSVs:
   i = 0;
   
   for (Index za_sca = 0; za_sca < scat_za_grid.nelem(); za_sca ++)
     {
       for (Index aa_sca = 0; aa_sca < scat_aa_grid.nelem(); aa_sca ++)
         {
           for (Index za_inc = 0; za_inc < scat_za_grid.nelem(); za_inc ++)
             {
               for (Index aa_inc = 0; aa_inc < scat_aa_grid.nelem(); aa_inc ++)
                 {
                   scat_theta_gps[za_sca][aa_sca][za_inc][aa_inc] = 
                     gps_array[i];
                   scat_theta_itws(za_sca, aa_sca, za_inc, aa_inc, joker) =
                     itws_array(i, joker);
                   i++;
                   //                   cout << "gps " <<  scat_theta_gps[za_sca][aa_sca][za_inc][aa_inc] << endl;
                 }
             }
         }
     }
}


//! Prepare single scattering data for a DOIT scattering calculation.
/*! 
  This function has can be used for scatttering calcualtions using the 
  DOIT method. It is a method optimized for randomly oriented 
  scattering media and it can only be used for such cases. 

  It has to be used in *scat_mono_agenda*. The phase matrix data is 
  interpolated on the actual frequency and on all scattering angles 
  following from all possible combinations in *scat_za_grid* and 
  *scat_aa_grid*.
  
  Output:
  \param scat_theta All scattering angles.
  \param pha_mat_sptDOITOpt Interpolated phase matrix data.
    Input:
  \param scat_za_grid Zenith angle grid for scattering calculations.
  \param scat_aa_grid Azimuth angle grid for scattering calculations.
  \param scat_data_raw Array of single scattering data for all hydrometeor
                      species.
  \param f_grid Frequency grid.
  \param f_index Frequency index. 
  \param atmosphere_dim Atmospheric dimension. 

  \author Claudia Emde
  \date   2003-08-25
*/
void ScatteringDataPrepareDOITOpt( //Output:
                                   Tensor4& scat_theta,
                                   Tensor6& pha_mat_sptDOITOpt,
                                   //Input:
                                   const Index& za_grid_size,
                                   const Vector& scat_aa_grid,
                                   const ArrayOfSingleScatteringData& scat_data_raw,
                                   const Vector& f_grid,
                                   const Index& f_index,
                                   const Index& atmosphere_dim
                                   )
{

  // For 1D calculation the scat_aa dimension is not required:
  Index N_aa_sca;
  if(atmosphere_dim == 1)
    N_aa_sca = 1;
  else
    N_aa_sca = scat_aa_grid.nelem();
  
  Vector za_grid;
  nlinspace(za_grid, 0, 180, za_grid_size);

  Index N_pt = scat_data_raw.nelem();  

  // Initialize variables:
  scat_theta.resize(za_grid.nelem(), N_aa_sca, 
                    za_grid.nelem(),scat_aa_grid.nelem());
  
  pha_mat_sptDOITOpt.resize(N_pt, za_grid_size, N_aa_sca, 
                    za_grid_size, scat_aa_grid.nelem(), 6);
  
  // Only one  particle type is taken.
  for( Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      // Frequency interpolation of the data:
      Matrix pha_mat_data_int(pha_mat_data_raw.nshelves(), 6);
      
      // Gridpositions:
      GridPos freq_gp;
      gridpos(freq_gp, f_datagrid, f_grid[f_index]); 
      
      // Interpolationweights:
      Vector freq_itw(2);
      interpweights(freq_itw, freq_gp);
  
  
      for (Index i = 0; i < 6; i++)
        {
          for (Index za_idx = 0; za_idx < za_datagrid.nelem(); za_idx ++)
            {
              pha_mat_data_int(za_idx, i) =
                interp(freq_itw, pha_mat_data_raw(joker, za_idx, 0, 0, 0, i),
                       freq_gp);
            }
        }
      
      // Calculate all scattering angles for all combinations of incoming 
      // and scattered directions and interpolation.
      
      Numeric za_sca_rad, za_inc_rad, aa_sca_rad, aa_inc_rad;
      Numeric theta_rad;
      Numeric theta;
      GridPos thet_gp;
      Vector itw(2);
      
      for (Index za_sca = 0; za_sca < za_grid_size; za_sca ++)
        {
          for (Index aa_sca = 0; aa_sca < N_aa_sca; aa_sca ++)
            {
              for (Index za_inc = 0; za_inc < za_grid_size; za_inc ++)
                {
                  for (Index aa_inc = 0; aa_inc < scat_aa_grid.nelem(); aa_inc ++)
                    {
                      za_sca_rad = za_grid[za_sca] * DEG2RAD;
                      za_inc_rad = za_grid[za_inc] * DEG2RAD;
                      aa_sca_rad = scat_aa_grid[aa_sca] * DEG2RAD;
                      aa_inc_rad = scat_aa_grid[aa_inc] * DEG2RAD;
                      
                      theta_rad = acos(cos(za_sca_rad) * cos(za_inc_rad) + 
                                       sin(za_sca_rad) * sin(za_inc_rad) * 
                                       cos(aa_sca_rad - aa_inc_rad));
                      
                      scat_theta(za_sca, aa_sca, za_inc, aa_inc)
                        = theta_rad;
                      
                      theta = theta_rad * RAD2DEG;
                      
                      gridpos(thet_gp, za_datagrid, theta);
                      interpweights(itw, thet_gp);
                      
                      for (Index i = 0; i < 6; i++)
                        {
                          pha_mat_sptDOITOpt(i_pt, za_sca, aa_sca, za_inc, aa_inc, i)
                            = interp(itw, pha_mat_data_int(joker, i), 
                                     thet_gp);
                        }
                    }
                }
            }
        }
    } 
 }

 
//! No preparation of single scattering data.
/*! 
  The parameters below are set to be empty. If this function is used 
  scattering angles, grid positions and interpolation weights are calculated
  inside the WSM *pha_mat_sptFromData*. 
  
  Output:
  \param scat_theta All scattering angles.
  \param scat_theta_gps Grid positions of scattering angles.
  \param scat_theta_itws Interpolation weights.
  
  \author Claudia Emde
  \date   2003-08-20
*/
void ScatteringDataPrepareOFF( //Output:
                               Tensor4& scat_theta,
                               ArrayOfArrayOfArrayOfArrayOfGridPos& 
                               scat_theta_gps,
                               Tensor5& scat_theta_itws
                               )
{ 
  scat_theta.resize(0,0,0,0);
  scat_theta_gps.resize(0);
  scat_theta_itws.resize(0,0,0,0,0);
}


//! scat_data_monoCalc
/*! 
  calculates monochromatic sngle scattering data for all particle types by
  interpolating scat_data_raw over frequency.

  \param scat_data_mono     Output: monochomatic scattering data
  \param scat_data_raw      raw scattering data WSV
  \param f_grid
  \param f_index

  \author Cory Davis
  \date 2003-12-16
*/
void scat_data_monoCalc(
			ArrayOfSingleScatteringData& scat_data_mono,
			const ArrayOfSingleScatteringData& scat_data_raw,
			const Vector& f_grid,
			const Index& f_index
			)
{
  const Index N_pt = scat_data_raw.nelem();
  
  //Initialise scat_data_mono
  scat_data_mono.resize(N_pt);
  
  // Loop over the included particle_types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {

      // Gridpositions:
      GridPos freq_gp;
      gridpos(freq_gp, f_datagrid, f_grid[f_index]); 

      // Interpolationweights:
      Vector itw(2);
      interpweights(itw, freq_gp);

      //Stuff that doesn't need interpolating
      scat_data_mono[i_pt].ptype=part_type;
      scat_data_mono[i_pt].f_grid.resize(1);
      scat_data_mono[i_pt].f_grid=f_grid[f_index];
      scat_data_mono[i_pt].T_grid=scat_data_raw[i_pt].T_grid;
      scat_data_mono[i_pt].za_grid=za_datagrid;
      scat_data_mono[i_pt].aa_grid=aa_datagrid;

      //Phase matrix data
      scat_data_mono[i_pt].pha_mat_data.resize(1,pha_mat_data_raw.nshelves(),
				  pha_mat_data_raw.nbooks(),pha_mat_data_raw.npages(),
				  pha_mat_data_raw.nrows(),pha_mat_data_raw.ncols());
     
      for (Index i_za_sca = 0; i_za_sca < pha_mat_data_raw.nshelves() ; i_za_sca++)
        {
          for (Index i_aa_sca = 0; i_aa_sca < pha_mat_data_raw.nbooks(); i_aa_sca++)
            {
              for (Index i_za_inc = 0; i_za_inc < pha_mat_data_raw.npages(); 
                   i_za_inc++)
                {
                  for (Index i_aa_inc = 0; i_aa_inc < pha_mat_data_raw.nrows(); 
                       i_aa_inc++)
                    {  
                      for (Index i = 0; i < pha_mat_data_raw.ncols(); i++)
                        {
                          scat_data_mono[i_pt].pha_mat_data(0,i_za_sca, 
                                                i_aa_sca, i_za_inc, 
                                                i_aa_inc, i) =
                            interp(itw,
                                   pha_mat_data_raw(joker,i_za_sca, 
                                                i_aa_sca, i_za_inc, 
                                                i_aa_inc, i),
                                   freq_gp);
                        }
                    }
                }
            }
        }
      //Extinction matrix data
      scat_data_mono[i_pt].ext_mat_data.resize(1,ext_mat_data_raw.npages(),
					       ext_mat_data_raw.nrows(),
					       ext_mat_data_raw.ncols());
      for (Index i_za_sca = 0; i_za_sca < ext_mat_data_raw.npages() ; i_za_sca++)
        {
          for(Index i_aa_sca = 0; i_aa_sca < ext_mat_data_raw.nrows(); i_aa_sca++)
            {
              //
              // Interpolation of extinction matrix:
              //
              for (Index i = 0; i < ext_mat_data_raw.ncols(); i++)
                {
                  scat_data_mono[i_pt].ext_mat_data(0,i_za_sca, i_aa_sca, i) =
                    interp(itw, ext_mat_data_raw(joker, i_za_sca, i_aa_sca, i),
                           freq_gp);
                }
            }
        }
      //Extinction matrix data
      scat_data_mono[i_pt].abs_vec_data.resize(1,abs_vec_data_raw.npages(),
					       abs_vec_data_raw.nrows(),
					       abs_vec_data_raw.ncols());
      for (Index i_za_sca = 0; i_za_sca < abs_vec_data_raw.npages() ; i_za_sca++)
        {
          for(Index i_aa_sca = 0; i_aa_sca < abs_vec_data_raw.nrows(); i_aa_sca++)
            {
              //
              // Interpolation of absorption vector:
              //
              for (Index i = 0; i < abs_vec_data_raw.ncols(); i++)
                {
                  scat_data_mono[i_pt].abs_vec_data(0,i_za_sca, i_aa_sca, i) =
                    interp(itw, abs_vec_data_raw(joker, i_za_sca, i_aa_sca, i),
                           freq_gp);
                }
            }
        }

    }
}



//! Calculates opticle properties for the single particle types.
/*! 
  This is essentially a monochromatic version of opt_prop_sptFromData. 
ext_mat_spt and abs_vec_spt are calculated from scat_data_mono instead of 
scat_data_raw

\param ext_mat_spt  Output: extinction matrix for each particle at the given
                         propagation direction
\param abs_vec_spt  Output: absorption coefficient vector at the given
                         propagation direction
\param scat_data_mono monochromatic single scattering data
\param scat_za_grid Workspace variable
\param scat_aa_grid Workspace variable
\param scat_za_index Workspace variable
\param scat_aa_index Workspace variable

\author Cory Davis
\date 2003-12-17

*/
void opt_prop_sptFromMonoData( // Output and Input:
                         Tensor3& ext_mat_spt,
                         Matrix& abs_vec_spt,
                         // Input:
                         const ArrayOfSingleScatteringData& scat_data_mono,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid,
                         const Index& scat_za_index, // propagation directions
                         const Index& scat_aa_index
                         )
{
  const Index N_pt = scat_data_mono.nelem();
  const Index stokes_dim = ext_mat_spt.ncols();
  const Numeric za_sca = scat_za_grid[scat_za_index];
  const Numeric aa_sca = scat_aa_grid[scat_aa_index];

  if (stokes_dim > 4 || stokes_dim < 1){
    throw runtime_error("The dimension of the stokes vector \n"
                         "must be 1,2,3 or 4");
  }
  
  assert( ext_mat_spt.npages() == N_pt );
  assert( abs_vec_spt.nrows() == N_pt );

  // Initialisation
  ext_mat_spt = 0.;
  abs_vec_spt = 0.;


  // Loop over the included particle_types
  for (Index i_pt = 0; i_pt < N_pt; i_pt++)
    {
      // First we have to transform the data from the coordinate system 
      // used in the database (depending on the kind of particle type 
      // specified by *ptype*) to the laboratory coordinate sytem. 
      
      //
      // Do the transformation into the laboratory coordinate system.
      //
      // Extinction matrix:
      //
     
  
      ext_matTransform(ext_mat_spt(i_pt, joker, joker),
                       scat_data_mono[i_pt].ext_mat_data(0,joker,joker,joker),
                       scat_data_mono[i_pt].za_grid, 
		       scat_data_mono[i_pt].aa_grid, 
		       scat_data_mono[i_pt].ptype,
                       za_sca, aa_sca);
      // 
      // Absorption vector:
      //
      abs_vecTransform(abs_vec_spt(i_pt, joker),
                       scat_data_mono[i_pt].abs_vec_data(0,joker,joker,joker),
                       scat_data_mono[i_pt].za_grid, 
		       scat_data_mono[i_pt].aa_grid, 
		       scat_data_mono[i_pt].ptype,
                       za_sca, aa_sca);                
    }


}
                          
