/*!
  \file   workspace.cc
  \brief  Definition of function wsv_data.

  This file contains the function define_wsv_data, which
  sets the WSV group names and the lookup data for the WSVs.
  You have to edit this function whenever you add a new
  workspace variable. 

  \author Stefan Buehler
  \date   2000-06-10
*/
#include "workspace.h"
#include "tokval_io.h"

// Some #defines to make the records better readable:
#define NAME(x) x
#define DESCRIPTION(x) x
#define GROUP(x) x

namespace global_data {
Array<WsvRecord> wsv_data;

std::map<String, Index> WsvMap;
}  // namespace global_data

using global_data::wsv_data;

void define_wsv_data() {
  //--------------------< Build the wsv data >--------------------
  // Initialize to empty, just in case.
  wsv_data.resize(0);

  /* Templace record entry:

  wsv_data.push_back
    (WsvRecord
     ( NAME( "workspace_variable_name" ),
       DESCRIPTION
       (
        "Brief description of the variable (1 line).\n"
        "\n"
        "Detailed description of the variable. Don't be too short here,\n"
        "this is the main place where your documentation should be. I\n"
        "really recommend to edit this in a text buffer, so that you can\n"
        "do some re-formatting until it looks nice. Only at the end put it\n"
        "in quotes and add the line breaks.\n"
        "\n"
        "Use blank lines to separate paragraphs.  There really should be a\n"
        "detailed descriptions of all component of your variable, if it\n"
        "has a complicated type. Also some detailed discussion of the\n"
        "dimensions if necessary. Also some detailed discussion of the\n"
        "members if your variable is a structure.\n"
        "\n"
        "Usage:      Set by user (or "Method output.")\n"
        "\n"
        "Units:      E.g., kg/m\n"
        "\n"
        "Dimensions: [ first dimension, second dimension, ... ]\n"
        "or\n"
        "Size:       [ .., nrows, ncols ]\n"
        "\n"
        "Members:    Here you would list the members if your\n"
        "            variable is a structure.\n"
        "\n"
        "Dimensions: [x, y]\n"
        "\n"
        "Unit: Which unit this variable uses\n"
        "\n"
        "Give the keywords above only if they apply, i.e., Members only\n"
        "for a structure, Units only for a physical variable.\n"
        "Use either Dimensions or Size, depending on what is most appropiate\n"
        "for the variable.\n"
        ),
      GROUP( "VariableType" )));

*/

  /*----------------------------------------------------------------------
    Let's put in the variables in alphabetical order. This gives a clear
    rule for where to place a new variable and this gives a nicer
    results when the methods are listed by "arts -w all".  No
    distinction is made between uppercase and lowercase letters. The
    sign "_" comes after all letters.
    Patrick Eriksson 2002-05-08
  ----------------------------------------------------------------------*/

  wsv_data.push_back(WsvRecord(
      NAME("aa_grid"),
      DESCRIPTION(
          "Azimuthal angle grid.\n"
          "\n"
          "The azimutal angle grid, on which the *cloudbox_field* is stored. \n"
          "This grid is used for RT calculations inside the cloudbox, \n"
          "therefore one has to define it if the cloudbox is activated by \n"
          "the flag *cloudbox_on*. Furthermore the zenith angle grid is also used"
          "for RT calculations of clear-sky *spectral_radiance_field*.\n"
          "The grid must be sorted in increasing order, with no repetitions.\n"
          "\n"
          "Usage:      Set by the user.\n"
          "\n"
          "Unit:       degrees \n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("aa_index"),
      DESCRIPTION(
          "Azimuth angle index for scattering calculations.\n"
          "\n"
          "This variable is used in methods used for computing scattering\n"
          "properties. \n"
          "It holds the information about the azimuth angles for which the \n"
          "scattering calculations are done.  The angles used for computing \n"
          "scattering properties of particles can be different from that used \n"
          "for radiative transfer calculation. \n"
          "\n"
          "Usage:    Method output.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_cia_data"),
      DESCRIPTION(
          "HITRAN Collision Induced Absorption (CIA) Data.\n"
          "\n"
          "This variable holds HITRAN CIA data (binary absorption\n"
          "cross-sections). The data itself is described in: Richard, C. et al.\n"
          "(2012), New section of the HITRAN database: Collision-induced\n"
          "absorption (CIA), J. Quant. Spectrosc. Radiat. Transfer, 113,\n"
          "1276-1285, doi:10.1016/j.jqsrt.2011.11.004.\n"
          "\n"
          "The binary absorption cross-sections have to be multiplied with the\n"
          "densities of both molecules to get absorption coefficients.\n"
          "\n"
          "Dimensions:\n"
          "\n"
          "The outer array dimension in the ArrayOfArrayOfCIARecord is the same\n"
          "as that of *abs_species*. There will be CIA data only for those\n"
          "species that contain a CIA tag, for all other species it will be\n"
          "empty. The inner array dimension corresponds to the number of CIA tags\n"
          "for this species (there could be for example both N2-N2 and N2-H2) in\n"
          "the same species.\n"
          "\n"
          "The CIA *abs_species* tags are described in *abs_speciesSet*.\n"
          "\n"
          "Each individual CIARecord holds the complete information from one\n"
          "HITRAN CIA file. For the given pair of molecules A HITRAN CIA data\n"
          "file can hold several datasets (data for different temperatures but\n"
          "fixed frequency range).\n"
          "\n"
          "Units:\n\n"
          " - Frequencies: Hz\n"
          " - Binary absorption cross-sections: m^5*molecule^-2\n"),
      GROUP("ArrayOfCIARecord")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_f_interp_order"),
      DESCRIPTION(
          "Frequency interpolation order for absorption lookup table.\n"
          "\n"
          "The interpolation order to use when interpolating the absorption\n"
          "lookup table in frequency. This is in particular needed for\n"
          "calculations with Doppler shift, so that absorption is interpolated to\n"
          "the shifted frequency grid. One is linear interpolation, two\n"
          "quadratic, and so on.\n"
          "\n"
          "As a special case, order 0 in this particular case means no\n"
          "interpolation. In that case f_grid must match exactly the grid inside\n"
          "the lookup table. This is the global default value.\n"),
      GROUP("Index"), Index{0}));
  

  wsv_data.push_back(WsvRecord(
      NAME("abs_hitran_relmat_data"),
      DESCRIPTION(
          "HITRAN line mixing data to compute the relaxation matrix.\n"
          "\n"
          "This variable holds HITRAN line mixing data\n"
          "as per J. Lamouroux, L. Realia, X. Thomas, et al., J.Q.S.R.T. 151 (2015), 88-96\n"
          "\n"
          "It is used for absorption bands with these population tags:\n"
          " - ByHITRANFullRelmat\n"
          " - ByHITRANRosenkranzRelmat\n"),
      GROUP("HitranRelaxationMatrixData")));
  
  wsv_data.push_back(WsvRecord(NAME("abs_lines"),
                               DESCRIPTION("A list of spectral line data.\n"),
                               GROUP("ArrayOfAbsorptionLines")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_lines_per_species"),
      DESCRIPTION(
          "A list of spectral line data for each tag.\n"
          "\n"
          "Dimensions: [*abs_species*.nelem()][Depends on how many bands there are in *abs_lines*]\n"),
      GROUP("ArrayOfArrayOfAbsorptionLines")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_lookup"),
      DESCRIPTION(
          "An absorption lookup table.\n"
          "\n"
          "It holds an absorption lookup table, as well as all information that\n"
          "is necessary to use the table to extract absorption. Extraction\n"
          "routines are implemented as member functions. \n"
          "\n"
          "It has quite a complicated structure. For details see the Arts User\n"
          "Guide section \"The gas absorption lookup table\" or the source code\n"
          "documentation in gas_abs_lookup.h.\n"),
      GROUP("GasAbsLookup")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_nls"),
      DESCRIPTION(
          "Nonlinear species for absorption lookup table generation.\n"
          "\n"
          "A list of absorption species that should be treated non-linearly.\n"
          "This means that the H2O VMR should be varied when calculating the\n"
          "lookup table for those species.\n"
          "\n"
          "A typical example is for this to containt the Rosenkranz full\n"
          "absorption model species for water vapor and oxygen \n"
          "([\"H2O-PWR98\", \"O2-PWR93\"]).\n"
          "\n"
          "See user guide and online documentation of *abs_lookupCalc*\n"
          "for more details and usage examples.\n"),
      GROUP("ArrayOfArrayOfSpeciesTag")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_nls_pert"),
      DESCRIPTION(
          "Fractional perturbations for the nonlinear species in the absorption\n"
          "lookup table.\n"
          "\n"
          "This is a vector of fractional perturbations that should contain 1\n"
          "(the unperturbed reference profile). A value of 0 may lead to error\n"
          "messages from some absorption routines, so a possible content for this\n"
          "variable is: [1e-24, 1, 2].\n"
          "(This is similar to *abs_t_pert*, but multiplicative, not additive.)\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_nls_interp_order"),
      DESCRIPTION(
          "The interpolation order to use when interpolating absorption between\n"
          "the H2O values given by *abs_nls_pert*.\n"
          "\n"
          "This is used by methods extracting absorption coefficients\n"
          "from the lookup table, and by methods setting up\n"
          "parameters for lookup table generation.\n"
          "\n"
          "Note that the number of points used in the interpolation scheme is\n"
          "interpolation order + 1 (e.g., two for first order interpolation).\n"),
      GROUP("Index"), Index{5}));

  wsv_data.push_back(WsvRecord(
      NAME("abs_p_interp_order"),
      DESCRIPTION(
          "The interpolation order to use when interpolating absorption\n"
          "between pressure levels.\n"
          "\n"
          "This is used by methods extracting absorption coefficients\n"
          "from the lookup table, and by methods\n"
          "setting up parameters for lookup table generation.\n"
          "\n"
          "Note that the number of points used in the interpolation scheme is\n"
          "interpolation order + 1 (e.g., two for first order interpolation).\n"),
      GROUP("Index"), Index{5}));

  wsv_data.push_back(WsvRecord(
      NAME("abs_t_pert"),
      DESCRIPTION(
          "Temperature perturbations for the absorption lookup table.\n"
          "\n"
          "This is a vector containing temperature perturbations (in Kelvin) that\n"
          "should be added to the reference temperature profile. (Similar to\n"
          "*abs_nls_pert*, but additive, not multiplicative.) Should normally\n"
          "contain 0, to include the reference profile itself. Example content:\n"
          "[-5, 0, 5].\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_t_interp_order"),
      DESCRIPTION(
          "The interpolation order to use when interpolating absorption between\n"
          "the temperature values given by *abs_t_pert*.\n"
          "\n"
          "This is used by methods\n"
          "extracting absorption coefficients from the lookup table, and by\n"
          "methods setting up parameters for lookup table generation.\n"
          "\n"
          "Note that the number of points used in the interpolation scheme is\n"
          "interpolation order + 1 (e.g., two for first order interpolation).\n"),
      GROUP("Index"), Index{7}));

  wsv_data.push_back(WsvRecord(
      NAME("abs_lookup_is_adapted"),
      DESCRIPTION(
          "Flag to indicate whether *abs_lookupAdapt* has already been\n"
          "called.\n"
          "\n"
          "Values:\n\n - 0=false\n - 1=true.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_p"),
      DESCRIPTION(
          "List of pressures to be used for the calculation of absorption\n"
          "coefficients.\n"
          "\n"
          "This can be copied from the global *p_grid*, but could also be\n"
          "different.\n"
          "\n"
          "Any absorption method should check that the length of this vector\n"
          "is the same as that of *abs_t*\n"
          "\n"
          "Dimension: [p_grid]\n"
          "\n"
          "Unit: Pa\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_species"),
      DESCRIPTION(
          "Tag groups for gas absorption.\n"
          "\n"
          "This is an array of arrays of SpeciesTag tag definitions. It defines the\n"
          "available tag groups for the calculation of scalar gas absorption\n"
          "coefficients.  See online documentation of method *abs_speciesSet* for\n"
          "more detailed information how tag groups work and some examples.\n"),
      GROUP("ArrayOfArrayOfSpeciesTag")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_t"),
      DESCRIPTION(
          "List of temperatures to be used for the calculation of absorption\n"
          "coefficients.\n"
          "\n"
          "In contrast to the global *t_field*, this is just a vector. Any\n"
          "absorption method should check that the length of this vector is the\n"
          "same as that of *abs_p*\n"
          "\n"
          "Dimension: [p_grid]\n"
          "\n"
          "Unit: K\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_nlte"),
      DESCRIPTION(
          "NLTE temperatures or ratios to be used for the calculation of\n"
          "absorption coefficients.\n"
          "\n"
          "In contrast to the global *nlte_field*, this is just a matrix. Any\n"
          "absorption method should check that the columns of this vector is the\n"
          "same as that of *abs_p*\n"
          "\n"
          "Dimension: [nltes, 1, 1, p_grid] or [ 0, 0, 0, 0 ]\n"
          "\n"
          "Unit: K\n"),
      GROUP("EnergyLevelMap")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_vec"),
      DESCRIPTION(
          "Total absorption vector.\n"
          "\n"
          "This variable contains the absorption coefficient vector which\n"
          "is used in the RTE calculation. It is the physical absorption which\n"
          "includes particle absorption for all considered scattering elements as\n"
          "well as gaseous absorption for all selected gaseous species.\n"
          "The vector is calculated by *opt_prop_bulkCalc*\n"
          "The dimension of the variable adapts to *stokes_dim*.\n"
          "\n"
          "See ARTS user guide (AUG) for further information. Use the index to find\n"
          "where this variable is discussed. The variable is listed as a subentry\n"
          "to \"workspace variables\".\n"
          "\n"
          "Usage:      Output of *opt_prop_bulkCalc* \n"
          "\n"
          "Unit:       m^2\n"  //FIXME: really m2? not 1/m?
          "\n"
          "Dimensions: [f_grid, stokes_dim]\n"),
      GROUP("StokesVector")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_vec_spt"),
      DESCRIPTION(
          "Absorption vectors of the scattering elements.\n"
          "\n"
          "This variable contains the elements of the absorption vector of the\n"
          "individual scattering elements. It is calculated in the agenda \n"
          "*spt_calc_agenda*.\n"
          "\n"
          "See ARTS user guide (AUG) for further information.\n"
          "\n"
          "Usage:      Input and Output of the method abs_vec_sptCalc\n"
          "\n"
          "Unit:        m^2\n"  //FIXME: really m2? not 1/m?
          "\n"
          "Dimensions: [number of scattering elements, stokes_dim]\n"),
      GROUP("ArrayOfStokesVector")));

  wsv_data.push_back(WsvRecord(
      NAME("abs_vmrs"),
      DESCRIPTION("The VMRs (unit of absolute number) on the abs_p grid.\n"
                  "\n"
                  "Dimensions: [tag_groups.nelem(), abs_p.nelem()]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("agenda_array_index"),
      DESCRIPTION(
          "Index of the current agenda in *ArrayOfAgenda*.\n"
          "\n"
          "This is set during the execution of an agenda from an *ArrayOfAgenda*.\n"
          "It indicates the index of the current agenda inside the array.\n"
          "\n"
          "Unit:  Integer value.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("antenna_dim"),
      DESCRIPTION(
          "The dimensionality of the antenna pattern (1-2).\n"
          "\n"
          "A dimensionality of 1 means that only the respons variation in the\n"
          "zenith direction is considered. The provided respons shall then be the\n"
          "integrated in the azimuth direction. For 2D, the respons of the\n"
          "antenna has both a zenith and azimuth variation.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  Integer value [1-2].\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("antenna_dlos"),
      DESCRIPTION(
          "The relative line-of-sight of each antenna pattern.\n"
          "\n"
          "This variable describes the line-of-sight of the individual antennae\n"
          "relative to *sensor_los*. If each measurement block corresponds to\n"
          "a single antenna pattern, the normal choice is to set the angle(s) of\n"
          "this variable to zero.\n"
          "\n"
          "The first column holds the relative zenith angle. This column is\n"
          "mandatory for all atmospheric dimensionalities. For 3D, there can\n"
          "also be a second column, giving relative azimuth angles. If this\n"
          "column is not present (for 3D) zero azimuth off-sets are assumed.\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  [ degrees, degrees ]\n"
          "\n"
          "Size:  [ number of antennae, 1 or 2 ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("antenna_response"),
      DESCRIPTION(
          "The antenna pattern/response.\n"
          "\n"
          "This WSV describes the antenna response as a function of polarisation\n"
          "(pol), frequencue (f), zenith angle (za) and azimuth angle (aa).\n"
          "\n"
          "Polarisation dimension: If this dimension has size 1, the data are\n"
          "applied for all polarisations of concern. The data are otherwise used\n"
          "in sequential order. This signifies that, in general, the first\n"
          "polarisation \"layer\" corresponds to the first stokes dimension\n"
          "etc. An exception is if a polarisation rotation has been applied. In\n"
          "any case, it is up to the user to ensure that polarisations are\n"
          "consistently defined.\n"
          "\n"
          "Frequency dimension: If this dimension has size 1, the data are\n"
          "applied for all frequencies of concern. The given frequency must be\n"
          "inside the frequency range of concern. A linear interpolation is\n"
          "otherwise applied.\n"
          "\n"
          "Zenith angle dimension: This dimension must always have a size >= 2\n"
          "The response outside covered grid range is treated as zero. If\n"
          "*antenna_dim* equals 1, the data should correspond to the response\n"
          "integrated in the azimuthal direction.\n"
          "\n"
          "Azimuth angle dimension: If *antenna_dim* equals 1, this dimension\n"
          "must have size 1. A size >= 2 is otherwise required. The response\n"
          "outside covered grid range is treated as zero.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit: The unit of the actual response is 1/sr. Properly normalised\n"
          "a 4pi integration shall shall give 1.\n"
          "\n"
          "Dimensions: \n"
          "\n"
          "- GriddedField4:\n"
          "\n"
          "  - ArrayOfString field_names[N_pol]\n"
          "  - Vector f_grid[N_f]\n"
          "  - Vector za_grid[N_za]\n"
          "  - Vector aa_grid[N_aa]\n"
          "  - Tensor4 data[N_pol][N_f][N_za][N_aa]\n"),
      GROUP("GriddedField4")));

  wsv_data.push_back(WsvRecord(
      NAME("atmosphere_dim"),
      DESCRIPTION(
          "The atmospheric dimensionality (1-3).\n"
          "\n"
          "This variable defines the complexity of the atmospheric structure.\n"
          "The dimensionality is given by an integer between 1 and 3, where 1\n"
          "means 1D etc. This is the master variable for the atmospheric\n"
          "dimensionality, variables which size changes with the dimensionality\n"
          "are checked to match this variable. \n"
          "\n"
          "Methods adapt automatically to this variable. That is, it should\n"
          "not be needed to change any methods if the dimensionality is\n"
          "changed. However, not all methods are working for higher dimensions.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit: Integer value.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("atmfields_checked"),
      DESCRIPTION(
          "OK-flag for atmospheric grids and (physical) fields.\n"
          "\n"
          "The variable flags that clear-sky part of the atmosphere is\n"
          "defined in formally correct way. Example on problems captured\n"
          "include that the size of an atmospheric fields does not match the\n"
          "length of the atmospheric grids, and physically incorrect data such\n"
          "as negative temperatures.\n"
          "\n"
          "Note that *z_field* is not covered by this variable, it is instead\n"
          "treated to be part of the geometrical considerations where the ok-flag\n"
          "is denoted as *atmgeom_checked*. The cloudbox is covered by\n"
          "*cloudbox_checked*.\n"
          "\n"
          "Shall be set by *atmfields_checkedCalc*. See that WSMs for treated\n"
          "WSVs. Only the value 1 is taken as OK.\n"),
      GROUP("Index"), Index{0}));

  wsv_data.push_back(WsvRecord(
      NAME("atmgeom_checked"),
      DESCRIPTION(
          "OK-flag for the geometry of the model atmosphere.\n"
          "\n"
          "The variable flags that reference ellipsoid, the surfae and *z_field*\n"
          "contain formally correct values. Includes for example, that *z_field*\n"
          "holds strictly increasing values at each geographical position.\n"
          "\n"
          "See also *atmfields_checked*.\n"
          "\n"
          "Shall be set by *atmgeom_checkedCalc*. Only the value 1 is taken\n"
          "as OK.\n"),
      GROUP("Index"), Index{0}));

  wsv_data.push_back(WsvRecord(
      NAME("atm_fields_compact"),
      DESCRIPTION(
          "A compact set of atmospheric fields on a common set of grids.\n"
          "\n"
          "Data is supposed to contain basic atmsopheric fields for a RT\n"
          "calculation, i.e., temperature, altitude, and gas VMRs. It can\n"
          "furthermore contain fields describing scattering species like mass\n"
          "content, mass flux, number density of diverse scattering species.\n"
          "\n"
          "VMR fields are unitless, scattering species fields are supposed to be\n"
          "in SI units (i.e. kg/m3 for mass contents, kg/m2/s for mass flux,\n"
          "1/m3 for number densities).\n"
          "\n"
          "The data are stored in a *GriddedField4*.\n"
          "\n"
          "The first field in the matrix (i.e., first matrix column) has to be\n"
          "atmospheric pressure. Apart from this, the order of the fields is\n"
          "free. Field content (apart from pressure) is identified by their\n"
          "given field name tag. Furthermore, absorption species (e.g. VMR)\n"
          "fields and scattering species fields are related to *abs_species*\n"
          "and *scat_species* entries, respectively, by their field name tags.\n"
          "The tags must exhibit the following structure:\n"
          "\n"
          "0) species identifier:\n"
          "   Fields, supposed to be sorted into *vmr_field*, must be headed the\n"
          "   tag 'abs_species'. Names of scattering species fields likewise must\n"
          "   be headed by the 'scat_species' tag. Temperature and altitude\n"
          "   fields do not hold any heading tag.\n"
          "1) species name:\n"
          "   The (core) name of the field: 'T' for temperature, 'z' for\n"
          "   altitude, the absorption species name (e.g. 'H2O, 'O3', etc.) for\n"
          "   absorption species, the scattering species name (e.g. 'IWC') for\n"
          "   scattering species. For scattering species, this part is matched\n"
          "   against the scattering species name part of the *scat_species*\n"
          "   tags.\n"
          "2) field type:\n"
          "   This has to be given for scattering species only, indicating the\n"
          "   type of the scattering species fields, i.e. 'mass_density',\n"
          "   'mass_flux', 'number_density', 'mean_mass'.\n"
          "\n"
          "Dashes ('-') serve as delimiter, separating the elements of each\n"
          "field name tag.\n"
          "\n"
          "Usage: Used inside batch calculations, to hold successive atmospheric\n"
          "states from an *ArrayOfGriddedField4*.\n"
          "\n"
          "Dimensions: \n"
          "\n"
          "- GriddedField4:\n"
          "\n"
          "  - ArrayOfString field_names[N_fields]\n"
          "  - Vector p_grid[N_p]\n"
          "  - Vector lat_grid[N_lat]\n"
          "  - Vector lon_grid[N_lon]\n"
          "  - Tensor4 data[N_fields][N_p][N_lat][N_lon]\n"),
      GROUP("GriddedField4")));

  wsv_data.push_back(WsvRecord(
      NAME("avk"),
      DESCRIPTION(
          "Averaging kernel matrix.\n"
          "\n"
          "This matrix is the partial derivative of the retrieved state vector\n"
          "with respect to the measurement vector (*y*).\n"
          "\n"
          "Usage: Used and set by inversion methods. \n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("backend_channel_response"),
      DESCRIPTION(
          "The response of each backend channel.\n"
          "\n"
          "The response is given as an *ArrayOfGriddedField1*. The grid consists of\n"
          "relative frequencies. These relative frequencies are added to \n"
          "*f_backend* to obtain the absolute frequency for each response value.\n"
          "The actual data are the response at each frequency grid point.\n"
          "\n"
          "There are here two options. If the array has length 1, the same\n"
          "response is applied for all channels. Accordingly, this assumes that\n"
          "all channels have the same response function. The second option is to\n"
          "specify the response for each channel seperately. This signifies that\n"
          "the *backend_channel_response* array has either 1 or n elements, where\n"
          "n is the length of *f_backend*\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Size:\n"
          "\n"
          "- Array[N_ch]\n"
          "\n"
          "  - GriddedField1\n"
          "\n"
          "    - [N_f]\n"
          "    - [N_f]\n"),
      GROUP("ArrayOfGriddedField1")));

  wsv_data.push_back(WsvRecord(
      NAME("backend_channel_response_multi"),
      DESCRIPTION(
          "As *backend_channel_response* but describes an instrument with\n"
          "muliple mixer/receiver chains.\n"
          "\n"
          "See *f_backend_multi* for when to use this variable and size\n"
          "constraints.\n"
          "\n"
          "Usage: Set by the user.\n "),
      GROUP("ArrayOfArrayOfGriddedField1")));

  wsv_data.push_back(WsvRecord(
      NAME("batch_atm_fields_compact"),
      DESCRIPTION(
          "An array of compact atmospheric states.\n"
          "\n"
          "This is used to hold a set of *atm_fields_compact* for batch\n"
          "calculations. For further information see *atm_fields_compact*.\n"),
      GROUP("ArrayOfGriddedField4")));

  wsv_data.push_back(WsvRecord(
      NAME("band_identifiers"),
      DESCRIPTION(
          "An array of identifiers for bands.\n"
          "\n"
          "Used by line mixing calculations to identify which bands to match to the\n"
          "line database.\n"),
      GROUP("ArrayOfQuantumIdentifier")));

  wsv_data.push_back(WsvRecord(
      NAME("batch_cloudbox_limits"),
      DESCRIPTION("An array of *cloudbox_limits*.\n"
                  "\n"
                  "This is used to hold a set of *cloudbox_limits* for batch\n"
                  "calculations. \n"),
      GROUP("ArrayOfArrayOfIndex")));

  wsv_data.push_back(WsvRecord(
      NAME("batch_pnd_fields"),
      DESCRIPTION("An array of compact pnd states.\n"
                  "\n"
                  "This is used to hold a set of 1D *pnd_field* for batch\n"
                  "calculations. \n"),
      GROUP("ArrayOfTensor4")));

  wsv_data.push_back(WsvRecord(
      NAME("channel2fgrid_indexes"),
      DESCRIPTION(
          "Definition of backend frequency response, link to *f_grid*.\n"
          "\n"
          "The WSV is used to describe the frequency response of backend channels\n"
          "together with the accompanying WSV *channel2fgrid_weights*.\n"
          "\n"
          "This WSV links each channel to the elements of *f_grid*. In short it\n"
          "lists what elements of *f_grid* that are relevant for each channel.\n"
          "\n"
          "More precisely, the first dimension gives the number of output channels.\n"
          "Each ArrayOfIndex gives the index of the values in *f_grid* associated\n"
          "with the channel of concern. For a pure double-sideband receiver, where\n"
          "there is one monochromatic frequency per passband, this argument could\n"
          "look like: [[0,5],[1,4],[2,3],[7,8],[7,8]].\n"),
      GROUP("ArrayOfArrayOfIndex")));

  wsv_data.push_back(WsvRecord(
      NAME("channel2fgrid_weights"),
      DESCRIPTION(
          "Definition of backend frequency response, weighting of *f_grid*.\n"
          "\n"
          "The WSV is used to describe the frequency response of backend channels\n"
          "together with the accompanying WSV *channel2fgrid_indexes*.\n"
          "\n"
          "This WSV shall have excatly the same sizes as *channel2fgrid_indexes*.\n"
          "Each element gives the weight to be assigned to the associated\n"
          "monochromatic frequency. \n"),
      GROUP("ArrayOfVector")));

  wsv_data.push_back(WsvRecord(
      NAME("cloudbox_checked"),
      DESCRIPTION(
          "OK-flag for variables associated with the cloudbox.\n"
          "\n"
          "This variable flags that cloudbox variables are defined in a formally\n"
          "and practically correct way. For example, that there is sufficient\n"
          "space between the cloudbox and edges of the model atmosphere (for\n"
          "2D and 3D). Pure clear-sky variables are covered by\n"
          "*atmfields_checked* (and *atmgeom_checked*).\n"
          "\n"
          "Relevant checks are performed by *cloudbox_checkedCalc*. Only the\n"
          "value 1 is taken as OK.\n"),
      GROUP("Index"), Index{0}));

  wsv_data.push_back(WsvRecord(
      NAME("cloudbox_field"),
      DESCRIPTION(
          "The spectral radiance field inside the cloudbx.\n"
          "\n"
          "This variable is used to store the radiance field inside the cloud\n"
          "box, probably determined by a scattering solver method.\n"
          "\n"
          "That is, this variable matches *spectral_radiance_field* but holds\n"
          "a field that is restricted to the cloud box.\n"
          "\n"
          "Unit: W / (m^2 Hz sr) for each Stokes component.\n"
          "\n"
          "Size:\n"
          " [ f_grid,\n"
          " p_grid, \n"
          " lat_grid, \n"
          " lon_grid, \n"
          " za_grid,\n"
          " aa_grid,\n"
          " stokes_dim ]\n"
          "\n"
          "Note: For 1D, the size of the latitude, longitude and azimuth\n"
          "dimension (N_aa) are all 1.\n"),
      GROUP("Tensor7")));

  wsv_data.push_back(WsvRecord(
      NAME("cloudbox_field_mono"),
      DESCRIPTION(
          "Monochromatic radiation field inside the cloudbox.\n"
          "\n"
          "This variable is used to store the monochromatic radiation field \n"
          "inside the cloudbox which is found by an iterative solution (DOIT).\n"
          "Refer to AUG for further information.\n"
          "\n"
          "Usage: Method output. \n"
          "\n"
          "Unit: W / (m^2 Hz sr) for each Stokes component.\n"
          "\n"
          "Size:\n"
          "       [(cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
          "       (cloudbox_limits[3] - cloudbox_limits[2]) +1, \n"
          "       (cloudbox_limits[5] - cloudbox_limits[4]) +1, \n"
          "       N_za, N_aa, N_i ]\n"
          "\n"
          "Note: For 1D, the size of the azimuth angle dimension (N_aa) is\n"
          "always 1.\n"),
      GROUP("Tensor6")));

  wsv_data.push_back(WsvRecord(
      NAME("cloudbox_field_mono_old"),
      DESCRIPTION(
          "As *cloudbox_field_mono* but from previous iteration.\n"
          "\n"
          "This variable is used to store the intensity field inside the\n"
          "cloudbox while performing the iteration. One has to store the\n"
          "intensity field of the previous iteration to be able to do the \n"
          "convergence test after each iteration.\n"
          "Refer to AUG for more information.\n"
          "\n"
          "Usage: Method output. \n"
          "\n"
          "Unit: W / (m^2 Hz sr) for each Stokes component.\n"
          "\n"
          "Size:\n"
          "      [(cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
          "      (cloudbox_limits[3] - cloudbox_limits[2]) +1, \n"
          "      (cloudbox_limits[5] - cloudbox_limits[4]) +1, \n"
          "      N_za, N_aa, N_i ]\n"),
      GROUP("Tensor6")));

  wsv_data.push_back(WsvRecord(
      NAME("cloudbox_limits"),
      DESCRIPTION(
          "The limits of the cloud box.\n"
          "\n"
          "This variable defines the extension of the cloud box. The cloud box \n"
          "is defined to be rectangular in the used coordinate system, with \n"
          "limits exactly at points of the involved grids. This means, for \n"
          "example, that the vertical limits of the cloud box are two pressure \n"
          "levels. For 2D, the angular extension of the cloud box is between \n"
          "two points of the latitude grid, and likewise for 3D but then also \n"
          "with a longitude extension between two grid points. The latitude and\n"
          "longitude limits for the cloud box cannot be placed at the end \n"
          "points of the corresponding grid as it must be possible to calculate\n"
          "the incoming intensity field.\n"
          "\n"
          "The variable *cloudbox_limits* is an array of index value with\n"
          "length twice *atmosphere_dim*. For each dimension there is a lower \n"
          "limit and an upper limit. The order of the dimensions is as usual \n"
          "pressure, latitude and longitude. The upper limit index must be \n"
          "greater then the lower limit index. For example, \n"
          "*cloudbox_limits* = [0 5 4 11 4 11] means that cloud box extends\n"
          "between pressure levels 0 and 5, and latitude and longitude points 4\n"
          "and 11.\n"
          "\n"
          "If *cloudbox_on* = 0, the content of this variable is neglected, but\n"
          "it must be initiated to some dummy values.\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage: Set by the user, either directly or using a method\n"
          "checking the extension of scattering particles.\n"
          "\n"
          "Unit:  Index values.\n"
          "\n"
          "Size:  [ 2 * atmosphere_dim ]\n"),
      GROUP("ArrayOfIndex")));

  wsv_data.push_back(WsvRecord(
      NAME("cloudbox_on"),
      DESCRIPTION(
          "Flag to activate the cloud box.\n"
          "\n"
          "Scattering calculations are confined to a part of the atmosphere\n"
          "denoted as the cloud box. The extension of the cloud box is given by\n"
          "*cloudbox_limits*. This variable tells methods if a cloud box is\n"
          "activated or not. \n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  Boolean.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("complex_refr_index"),
      DESCRIPTION(
          "Complex refractive index (n) data.\n"
          "\n"
          "The variable works as a lookup-table of complex refractive index.\n"
          "The matter type (water, ice ...) is unspecified, it is up to the\n"
          "user to fill the variable with data for the expected matter.\n"
          "This variable type can be used to describe n of both the surface and\n"
          "atmospheric particles. For the surface, a dedicated variable exists:\n"
          "*surface_complex_refr_index*.\n"
          "\n"
          "The column dimension has always size 2, where the first and second\n"
          "column holds the real and imaginary part of n, respectively. The row\n"
          "dimension matches temperature, and the page dimension is frequency.\n"
          "Both the temperature and frequency dimensions grids are allowed to\n"
          "have length 1, which is interpreted as n being constant in that\n"
          "dimension.\n"
          "\n"
          "When mapping these data to the required frequencies and temperatures\n"
          "a bi-linear interpolation is applied.\n"
          "\n"
          "Unit:       -\n"
          "\n"
          "Dimensions: \n"
          " - Vector f_grid[N_f]\n"
          " - Vector T_grid[N_T]\n"
          " - ArrayOfString Complex[2]\n"
          " - Tensor3 data[N_f][N_T][2]\n"),
      GROUP("GriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("covmat_block"),
      DESCRIPTION(
          "Holds matrices used to set blocks in *covmat_sx* and *covmat_se*.\n"
          "\n"
          "The matrix contained in this block will be added to the blocks in\n"
          "in *covmat_sx* and *covmat_se* by the corresponding WSMs. Its dimensions\n"
          "must agree with gridpoints of the correlated retrieval quantities.\n"
          "\n"
          "Usage:   Used by the retrievalAdd functions.\n"),
      GROUP("Sparse")));

  wsv_data.push_back(WsvRecord(
      NAME("covmat_inv_block"),
      DESCRIPTION(
          "Holds matrices used to set the inverse blocks in *covmat_sx* and *covmat_se*.\n"
          "\n"
          "The matrix contained in this block will be used as the inverse of the matrix\n"
          "contained in covmat_block.\n"
          "\n"
          "Usage:   Used by the retrievalAdd functions.\n"),
      GROUP("Sparse")));

  wsv_data.push_back(WsvRecord(
      NAME("covmat_se"),
      DESCRIPTION(
          "Covariance matrix for observation uncertainties.\n"
          "\n"
          "This matrix (Se) describes the uncertainty of the measurement vector (*y*),\n"
          "and can be writtenn as::\n"
          "\n"
          "  Se = Seps + Kb * Sb * Kb'\n"
          "\n"
          "where Seps describes direct measurement errors (such as thermal noise),\n"
          "Kb is Jacobian for forward model parameters, and Sb describes the uncertainty\n"
          "of the forwatrd model parameters.\n"
          "\n"
          "Usage:   Used by inversion methods.\n"
          "\n"
          "Dimensions: [ y, y ]\n"),
      GROUP("CovarianceMatrix")));

  wsv_data.push_back(WsvRecord(
      NAME("covmat_sx"),
      DESCRIPTION(
          "Covariance matrix of a priori distribution\n"
          "\n"
          "This covariance matrix describes the Gaussian a priori distribution\n"
          "for an OEM retrieval. It is represented using a symmetric block matrix.\n"
          "covmat_sx can be used in two ways: Either with a block for each retrieval\n"
          "quantity or with a single block containing the full covariance matrix.\n"
          "\n"
          "Using a single block for each retrieval quantity has is advantageous for\n"
          "if the retrieval quantities are assumed to be independent. In this case,\n"
          "the covariance blocks can be added separately for each quantity and will\n"
          "allow optimizing matrix multiplications and inverses required for the OEM\n"
          "calculation.\n"
          "\n"
          "The other case of using a single-block covariance matrix is supported\n"
          "for convenience as well.\n"
          "\n"
          "Usage:   Used by inversion methods.\n"
          "\n"
          "Dimensions: \n"
          "     [ x, x ]\n"),
      GROUP("CovarianceMatrix")));

  wsv_data.push_back(WsvRecord(
      NAME("covmat_so"),
      DESCRIPTION(
          "Covariance matrix describing the retrieval error due to uncertainties of\n"
          "the observation system.\n"
          "\n"
          "That is: So = G * Se * G', where G is the gain matrix (*dxdy*).\n"
          "\n"
          "Usage: Set by the covmat_soCalc workspace method to characterize the error.\n"
          "of a successful OEM calculation.\n"
          "\n"
          "Dimensions:\n"
          "    [x,x]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("covmat_ss"),
      DESCRIPTION(
          "Covariance matrix describing the retrieval error due to smoothing.\n"
          "\n"
          "That is: Ss = (A-I) * Sx * (A-I)', where A is the averaging kernel "
          "matrix (*avk*).\n"
          "\n"
          "Usage: Set by the covmat_ssCalc workspace method to characterize the.\n"
          "errors of a successful OEM calculation.\n"
          "\n"
          "Dimensions:\n"
          "    [x,x]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("depolarization_factor"),
      DESCRIPTION(
          "Depolarization factor for the scattered gas.\n"
          "\n"
          "The variable accounts for the anisotropy of the scatterer.\n"
          "It is the ratio of intensities parallel and perpendicular \n"
          "to the plan of scattering. A table of measured values is \n"
          "given by Penndorf (1957). Some values are: H2=0.02, N2=0.03\n"
          "O2=0.06, CO2=0.09 and atmospheric air=0.03."),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("disort_aux"),
      DESCRIPTION(
          "Auxilary data to the output of the cloudbox_fieldDisort-Methods.\n"
          "\n"
          "Different data beside the direct result of Disort\n"
          "calculations can be obtained by this variable. These auxilary\n"
          "data are selected by *disort_aux_vars*.\n"
          "\n"
          "Usage:      Provided by some radiative transfer methods.\n"
          "\n"
          "Dimensions: [quantity][ f_grid, number of disort levels/layers ]\n"),
      GROUP("ArrayOfMatrix")));

  wsv_data.push_back(WsvRecord(
      NAME("disort_aux_vars"),
      DESCRIPTION(
          "Selection of quantities for *disort_aux*.\n"
          "\n"
          "Each element of this string array determines the quantity for the\n"
          "corresponding element in *disort_aux* (i.e. the quantities\n"
          "are stored in the order given in *disort_aux_vars*).\n"
          "\n"
          "The possible choices vary between the Disort methods. See the WSM you select\n"),
      GROUP("ArrayOfString"), ArrayOfString{}));


  wsv_data.push_back(WsvRecord(
      NAME("dlos"),
      DESCRIPTION(
          "A set of relative angles.\n"
          "\n"
          "This variable is a matrix having two columns. The two columns hold\n"
          "relative zenith angle and relative azimuth angle, respectively.\n"
          "\n"
          "These relative angles have zenith angle 90 deg as reference. This\n"
          "means that dza and daa represent the same angle distance at dza=0.\n"
          "\n"
          "Let us say that you add relative angles to a line-of-sight of za = 90\n"
          "and aa=0. Then adding the following (dza,daa) gives los of (za,aa):\n"
          "\n"
          "- (1,0) -> (91,0)\n"
          "- (0,1) -> (90,1)\n"
          "- (-90,45) -> (0,undefined)\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("dlos_weight_vector"),
      DESCRIPTION(
          "A weight associated with each direction *dlos*.\n"
          "\n"
          "A standard application should be to store the solid angle each\n"
          "row in *dlos* covers.\n"),
      GROUP("Vector")));
  
  wsv_data.push_back(WsvRecord(
      NAME("dobatch_calc_agenda"),
      DESCRIPTION(
          "Agenda defining the calculations to perform for each batch case.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("dobatch_cloudbox_field"),
      DESCRIPTION(
          "Batch of radiation fields.\n"
          "\n"
          "Each element of *dobatch_cloudbox_field* corresponds to a radiation field.\n"
          "See further *DOBatchCalc*.\n"
          "\n"
          "Usage: Most commonly produced by *DOBatchCalc*.\n"
          "\n"
          "Unit:  See *cloudbox_field*.\n"
          "\n"
          "Dimensions: Number of array elements equals number of batch cases.\n"),
      GROUP("ArrayOfTensor7")));

  wsv_data.push_back(WsvRecord(
      NAME("dobatch_radiance_field"),
      DESCRIPTION(
          "Batch of radiance fields.\n"
          "\n"
          "Each element of *dobatch_radiance_field* corresponds to a radiance field.\n"
          "See further *DOBatchCalc*.\n"
          "\n"
          "Usage: Most commonly produced by *DOBatchCalc*.\n"
          "\n"
          "Unit:  See *radiance_field*.\n"
          "\n"
          "Dimensions: Number of array elements equals number of batch cases.\n"),
      GROUP("ArrayOfTensor5")));

  wsv_data.push_back(WsvRecord(
      NAME("dobatch_irradiance_field"),
      DESCRIPTION(
          "Batch of irradiance fields.\n"
          "\n"
          "Each element of *dobatch_irradiance_field* corresponds to a irradiance field.\n"
          "See further *DOBatchCalc*.\n"
          "\n"
          "Usage: Most commonly produced by *DOBatchCalc*.\n"
          "\n"
          "Unit:  See *irradiance_field*.\n"
          "\n"
          "Dimensions: Number of array elements equals number of batch cases.\n"),
      GROUP("ArrayOfTensor4")));

  wsv_data.push_back(WsvRecord(
      NAME("dobatch_spectral_irradiance_field"),
      DESCRIPTION(
          "Batch of spectral irradiance fields.\n"
          "\n"
          "Each element of *dobatch_spectral_irradiance_field* corresponds to a\n"
          "spectral irradiance field.\n"
          "See further *DOBatchCalc*.\n"
          "\n"
          "Usage: Most commonly produced by *DOBatchCalc*.\n"
          "\n"
          "Unit:  See *spectral_irradiance_field*.\n"
          "\n"
          "Dimensions: Number of array elements equals number of batch cases.\n"),
      GROUP("ArrayOfTensor5")));

  wsv_data.push_back(WsvRecord(
      NAME("diy_dx"),
      DESCRIPTION(
          "Derivative of *iy* with respect to retrieval quantities.\n"
          "\n"
          "The variable gives the derivative if *iy* with respect to some\n"
          "variables (but not all jacobian variables). Handled are only variables\n"
          "affecting monochromatic pencil beam radiances where an (semi-)\n"
          "analytical expression can be applied (and that this calculation way\n"
          "has been selected when the jacobian was been set-up).\n"
          "\n"
          "The values in *diy_dx* considers the retrieval unit selected (such as\n"
          "\"nd\"), but no transformations are applied.\n"
          "\n"
          "Usage:      Output of *iy_main_agenda*.\n"
          "\n"
          "Dimensions: \n"
          "     [n_quantities][ n_retrieval_points, f_grid, stokes_dim ]\n"),
      GROUP("ArrayOfTensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("dpnd_data_dx"),
      DESCRIPTION(
          "Partial derivates of *pnd_data*.\n"
          "\n"
          "The variable gives the particle derivate of *pnd_data* with respect\n"
          "to the quantities set in *dpnd_data_dx_names*.\n"
          "\n"
          "Dimensions: [ n_quantities, n_points, n_scattering_elements ]\n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("dpnd_data_dx_names"),
      DESCRIPTION(
          "Selection of partial derivatives of *pnd_data*.\n"
          "\n"
          "This variable tells an element in *pnd_agenda_array* for which\n"
          "quantities partial derivatives shall be calculated.\n"
          "\n"
          "Dimensions: [ n_quantities ]\n"),
      GROUP("ArrayOfString")));

  wsv_data.push_back(WsvRecord(
      NAME("dpnd_field_dx"),
      DESCRIPTION(
          "Partial derivatives of *pnd_field*.\n"
          "\n"
          "The variable gives the particle derivative of *pnd_field* with respect\n"
          "to scattering species variables included in *jacobian_quantities*.\n"
          "\n"
          "The length of this array shall match the size of *jacobian_quantities*.\n"
          "For retrieval quantities that are not scattering species, the matching\n"
          "Tensor4 is of no relevance and must be set to be empty.\n"
          "\n"
          "Dimensions: [n_quantities][ n_scattering_elements, n_p, n_lat, n_lon ]\n"),
      GROUP("ArrayOfTensor4"), ArrayOfTensor4{}));

  wsv_data.push_back(WsvRecord(
      NAME("dpropmat_clearsky_dx"),
      DESCRIPTION(
          // FIXMEDOC
          "Partial derivative of absorption coefficients.\n"
          "\n"
          "This contains the partial derivative of absorption coefficients for\n"
          "one point in the atmosphere (one set of pressure, temperature, zn"
          "magnetic field, and VMR values) with respect to one of the input\n"
          "parameters.\n"
          "\n"
          "Dimension: [ n_quantities ] [naa, nza, nf, f(stokes_dim)]\n"
          "\n"
          "*jacobian_quantities* should be used to set the input variable for\n"
          "partial derivation\n"
          "\n"
          "Unit: 1/m/jacobian_quantity\n"),
      GROUP("ArrayOfPropagationMatrix")));

  wsv_data.push_back(WsvRecord(
      NAME("dpsd_data_dx"),
      DESCRIPTION(
          "Partial derivates of *psd_data*.\n"
          "\n"
          "The variable gives the particle derivate of *psd_data* with respect\n"
          "to the quantities set in *dpnd_data_dx_names*.\n"
          "\n"
          "Dimensions: [ n_quantities, n_points, n_scattering_elements ]\n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("dnlte_source_dx"),
      DESCRIPTION(
          "NLTE partial derivatives output is two parts:  S * dB/dx + dS/dx * B.\n"
          "\n"
          "Dimensions: [ quantities ] [nza, naa, nf, stokes_dim] or [0]\n"
          "\n"
          "Unit: 1/m/jacobian_quantity\n"),
      GROUP("ArrayOfStokesVector")));

  wsv_data.push_back(WsvRecord(
      NAME("doit_conv_flag"),
      DESCRIPTION(
          "Flag for the convergence test.\n"
          "\n"
          "This variable is initialized with 0 inside the method \n"
          "*cloudbox_field_monoIterate*.\n"
          "If after an iteration the convergence test is fulfilled, 1 is \n"
          "assigned which means that the iteration is completed. \n"
          "\n"
          "Usage: Method output. \n"),
      GROUP("Index")));

  wsv_data.push_back(
      WsvRecord(NAME("doit_conv_test_agenda"),
                DESCRIPTION("Agenda executing the DOIT convergence test.\n"),
                GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("doit_is_initialized"),
      DESCRIPTION("Flag to determine if *DoitInit* was called.\n"
                  "\n"
                  "This flag is checked by *DoitCalc* to make sure that\n"
                  "*DoitInit* was called before.\n"),
      GROUP("Index")));

  wsv_data.push_back(
      WsvRecord(NAME("doit_iteration_counter"),
                DESCRIPTION("Counter for number of iterations.\n"
                            "\n"
                            "This variable holds the number of iterations \n"
                            "while solving the VRTE using the DOIT method. \n"),
                GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("doit_mono_agenda"),
      DESCRIPTION("Agenda performing monochromatic DOIT calculation.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("doit_rte_agenda"),
      DESCRIPTION(
          "Agenda performing the DOIT cloudbox radiative transfer update.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("doit_scat_field_agenda"),
      DESCRIPTION(
          "Agenda calculating the scattering integral field in DOIT.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("doit_scat_field"),
      DESCRIPTION(
          "Scattered field inside the cloudbox.\n"
          "\n"
          "This variable holds the value of the scattering integral for all\n"
          "points inside the cloudbox. For more information refer to AUG.\n"
          "\n"
          "Usage: Input to ``cloudbox_fieldUpdate...``. \n"
          "\n"
          "Unit: W / (m^2 Hz sr) for each Stokes component.\n"
          "\n"
          "Size:\n"
          "     [(cloudbox_limits[1] - cloudbox_limits[0]) +1, \n"
          "     (cloudbox_limits[3] - cloudbox_limits[2]) +1, \n"
          "     (cloudbox_limits[5] - cloudbox_limits[4]) +1, \n"
          "     N_za, N_aa, N_i ]\n"),
      GROUP("Tensor6")));

  wsv_data.push_back(
      WsvRecord(NAME("doit_za_grid_opt"),
                DESCRIPTION("Optimized zenith angle grid.\n"
                            "\n"
                            "Output of the method *doit_za_grid_optCalc*.\n"
                            "\n"
                            "Usage:   Output of *doit_za_grid_optCalc*\n"
                            "\n"
                            "Unit:    degrees \n"),
                GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("doit_za_grid_size"),
      DESCRIPTION(
          "Number of equidistant grid points of the zenith angle grid.\n"
          "\n"
          "Grid points are defined from 0 to 180 deg, for the scattering\n"
          "integral calculation.\n"
          "\n"
          "Usage: Output of *DOAngularGridsSet*.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("doit_za_interp"),
      DESCRIPTION("Flag for interplation method in zenith angle dimension.\n"
                  "\n"
                  " - 0 - linear interpolation \n"
                  " - 1 - cubic interpolation \n"
                  "\n"
                  "Usage: Set by user in *doit_za_interpSet*. \n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("dsurface_emission_dx"),
      DESCRIPTION(
          "The derivative of *surface_emission* with respect to quantities\n"
          "listed in *dsurface_names*.\n"
          "\n"
          "Usage: Used internally of radiative transfer methods\n"
          "\n"
          "Dimensions: [dsurface_names][f_grid, stokes_dim]\n"),
      GROUP("ArrayOfMatrix")));

  wsv_data.push_back(WsvRecord(
      NAME("dsurface_names"),
      DESCRIPTION("Name of surface retrieval quantities.\n"
                  "\n"
                  "Usage: Used internally of radiative transfer methods\n"
                  "\n"
                  "Dimensions: [retrieval quantity]\n"),
      GROUP("ArrayOfString")));

  wsv_data.push_back(WsvRecord(
      NAME("dsurface_rmatrix_dx"),
      DESCRIPTION(
          "The derivative of *surface_rmatrix* with respect to quantities\n"
          "listed in *dsurface_names*.\n"
          "\n"
          "Usage: Used internally of radiative transfer methods\n"
          "\n"
          "Dimensions: [dsurface_names][surface_los, f_grid, stokes_dim, stokes_dim]\n"),
      GROUP("ArrayOfTensor4")));

  wsv_data.push_back(WsvRecord(
      NAME("dxdy"),
      DESCRIPTION(
          "Contribution function (or gain) matrix.\n"
          "\n"
          "This matrix is the partial derivative of the retrieved state vector\n"
          "with respect to the measurement vector (*y*).\n"
          "\n"
          "Usage: Used and set by inversion methods. \n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("ecs_data"),
      DESCRIPTION(
          "Error corrected sudden data\n"
          "\n"
          "Dimensions: [num IDs] [num Species]\n"
          "\n"
          "It is used for absorption bands with these population tags:\n"
          " - ByMakarovFullRelmat \n"
          " - ByRovibLinearDipoleLineMixing \n"),
      GROUP("MapOfErrorCorrectedSuddenData")));

  wsv_data.push_back(WsvRecord(
      NAME("ext_mat"),
      DESCRIPTION(
          "Total extinction matrix.\n"
          "\n"
          "This variable contains the extinction coefficient matrix, which\n"
          "is used in the RT calculation in the cloudbox. It is the physical\n"
          "extinction matrix which includes particle extinction for all chosen\n"
          "scattering species and gaseous extinction for all chosen gaseous species.\n"
          "\n"
          "See the ARTS user guide (AUG) for further information. Use the index to\n"
          "find where this variable is discussed. The variable is listed as a\n"
          "subentry to \"workspace variables\".\n"
          "\n"
          "Usage:      Output of *opt_prop_bulkCalc* \n"
          "\n"
          "Unit:       m^2\n"  //FIXME: really m2? not 1/m?
          "\n"
          "Dimensions: [f_grid, stokes_dim, stokes_dim]\n"),
      GROUP("PropagationMatrix")));

  wsv_data.push_back(WsvRecord(
      NAME("ext_mat_spt"),
      DESCRIPTION(
          "Extinction matrix for all individual scattering elements.\n"
          "\n"
          "This variable contains the elements of the extinction matrix of all\n"
          "individual scattering elements for a given propagation direction. It is\n"
          "calculated input as well as the output of the agenda *spt_calc_agenda*.\n"
          "\n"
          "Usage:      Output of *spt_calc_agenda* \n"
          "\n"
          "Unit:       m^2\n"  //FIXME: really m2? not 1/m?
          "\n"
          "Dimensions: [number of scattering elements, stokes_dim, stokes_dim]\n"),
      GROUP("ArrayOfPropagationMatrix")));

  wsv_data.push_back(WsvRecord(
      NAME("file_index"),
      DESCRIPTION(
          "Index number for files.\n"
          "\n"
          "See *WriteXMLIndexed* for further information.\n"
          "\n"
          "Usage:   Input to *WriteXMLIndexed* and *ReadXMLIndexed*. \n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(NAME("forloop_agenda"),
                               DESCRIPTION("Agenda performing a for loop.\n"),
                               GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("forloop_index"),
      DESCRIPTION(
          "The index for for-loops.\n"
          "\n"
          "This is the index that is used by method *ForLoop* to loop over\n"
          "*forloop_agenda*. \n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(NAME("fos_iyin_za_angles"),
                               DESCRIPTION("So far just testing of FOS ...\n"),
                               GROUP("Vector")));

  wsv_data.push_back(WsvRecord(NAME("fos_scatint_angles"),
                               DESCRIPTION("So far just testing of FOS ...\n"),
                               GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("f_backend"),
      DESCRIPTION(
          "The frequency position of each backend (spectrometer) channel.\n"
          "\n"
          "Usage: Set by the user.\n "
          "\n"
          "Unit:  Hz\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("f_backend_multi"),
      DESCRIPTION(
          "As *f_backend* but describes an instrument with muliple\n"
          "mixer/receiver chains.\n"
          "\n"
          "This variable is needed when e.g. the receiver has several mixers\n"
          "or the the receiver measures several polarisation and the channels\n"
          "differ in position or response function. \n"
          "\n"
          "The array has one element for each \"receiver chain\". The array\n"
          "length must match *backend_channel_response_multi*, and possibly\n"
          "also *lo_multi*.\n"
          "\n"
          "Usage: Set by the user.\n "
          "\n"
          "Unit:  Hz\n"),
      GROUP("ArrayOfVector")));

  wsv_data.push_back(WsvRecord(
      NAME("f_grid"),
      DESCRIPTION(
          "The frequency grid for monochromatic pencil beam calculations.\n"
          "\n"
          "Usage: Set by the user.\n "
          "\n"
          "Unit:  Hz\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("f_index"),
      DESCRIPTION(
          "Frequency index.\n"
          "\n"
          "Not all methods handle all monochromatic frequencies (of *f_grid*) in\n"
          "parellel and this variable is used for communication between methods,\n"
          "holding the index of the frequency treated presently.\n"
          "\n"
          "In some contexts, a negative f_index means all frequencies.\n"
          "\n"
          "Usage: Method output.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("gas_scattering_do"),
      DESCRIPTION(
          "Flag to activate gas scattering.\n"
          "\n"
          "If this variable is set to 0, no gas scattering will be considered,\n"
          "even if the gas_scattering_agenda is set.\n"
          "\n"),
      GROUP("Index"), Index{0}));

  wsv_data.push_back(WsvRecord(
      NAME("gas_scattering_output_type"),
      DESCRIPTION(
          "Flag to select the output of the *gas_scattering_agenda*.\n"
          "\n"
          "Internal communications variable, not intended to be used by user."
          "If equals 0 *gas_scattering_mat* is output and *gas_scattering_fct_legendre* is empty.\n"
          "If equals 1 *gas_scattering_fct_legendre* is output and *gas_scattering_mat* is empty.\n"
          "\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("gas_scattering_agenda"),
      DESCRIPTION(
          "Agenda calculating gas scattering extinction and phase matrix.\n"
          "\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("gas_scattering_los_in"),
      DESCRIPTION(
          "Incoming line-of-sight for gas scattering.\n"
          "\n"
          "This variable holds a local line-of-sight. The angles of this\n"
          "vector are defined as for *rte_los*.\n"
          "\n"
          "The WSV is used as input in *gas_scattering_agenda*\n"
          "\n"
          "Usage: Communication variable.\n"
          "\n"
          "Units: [ degree, degree ]\n"
          "\n"
          "Size:  [ 2 ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("gas_scattering_los_out"),
      DESCRIPTION(
          "Outgoing line-of-sight for gas scattering.\n"
          "\n"
          "This variable holds a local line-of-sight. The angles of this\n"
          "vector are defined as for *rte_los*.\n"
          "\n"
          "The WSV is used as input in *gas_scattering_agenda*\n"
          "\n"
          "Usage: Communication variable.\n"
          "\n"
          "Units: [ degree, degree ]\n"
          "\n"
          "Size:  [ 2 ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("gas_scattering_coef"),
      DESCRIPTION(
          "Spectrum of scattering coefficient matrices.\n"
          "\n"
          "This variable contains the elements of the extinction matrix solely\n"
          "due to scattering.\n"
          "\n"
          "Usage: Output of *gas_scattering_agenda*.\n"
          "\n"
          "Units: [ m^-1. ]\n"
          "\n"
          "Size:  [fgrid, stokes_dim, stokes_dim]\n"),
      GROUP("PropagationMatrix")));

  wsv_data.push_back(WsvRecord(
      NAME("gas_scattering_mat"),
      DESCRIPTION(
          "Spectrum of normalized phase matrices.\n"
          "\n"
          "This variable contains the elements of the normalized phase matrix\n"
          "for a specific incoming and outgoing direction.\n"
          "\n"
          "Usage: Output of *gas_scattering_agenda*.\n"
          "\n"
          "Units: [ 1 ]\n"
          "\n"
          "Size:  [fgrid, stokes_dim, stokes_dim]\n"),
      GROUP("TransmissionMatrix")));

  wsv_data.push_back(WsvRecord(
      NAME("gas_scattering_fct_legendre"),
      DESCRIPTION(
          "Normalized phase function as Legendre series.\n"
          "\n"
          "This variable contains the normalized phase function\n"
          "as Legendre series.\n"
          "\n"
          "Usage: Output of *gas_scattering_agenda*.\n"
          "\n"
          "Units: [ 1 ]\n"
          "\n"
          "Size:  [Number of Legendre polynomials]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("geo_pos"),
      DESCRIPTION(
          "Geo-position of a measurement.\n"
          "\n"
          "An empty vector is allowed, then flagging that no geo-positioning\n"
          "has been performed.\n"
          "\n"
          "Otherwise, this should be a vector having length 5. The elements are:\n"
          " - altitude\n"
          " - latitude\n"
          " - longitide\n"
          " - zenith angle\n"
          " - azimuth angle\n"
          "\n"
          "Dimensions: 0 or 5\n"
          "\n"
          "Unit:  [ m, deg, deg, deg, deg ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("g0"),
      DESCRIPTION(
          "Gravity at zero altitude.\n"
          "\n"
          "This variable is \"little g\" at the reference ellipsiod. That is,\n"
          "for Earth this is a value around 9.81 m/s2\n"),
      GROUP("Numeric")));

  wsv_data.push_back(
      WsvRecord(NAME("g0_agenda"),
                DESCRIPTION("Agenda providing the gravity constant.\n"),
                GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("heating_rates"),
      DESCRIPTION(
          "The heating rates of atmospheric layers.\n"
          "\n"
          "The heating rate is defined as the rate of temperature change of an \n"
          "atmospheric layer due to heating by absorption of radiation or if it\n"
          "is negative due to loss of energy by emission of radiation.\n"
          "\n"
          "Units: K s^-1\n"
          "\n"
          "Size: [ p_grid, lat_grid, lon_grid ]\n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("xsec_fit_data"),
      DESCRIPTION(
          "Fitting model coefficients for cross section species.\n"
          "\n"
          "Dimensions: [ n_species ]\n"
          "\n"
          "XsecRecord:\n"
          "\n"
          "- species: Name of species\n"
          "- version: Fit model version\n"
          "- fitcoeffs:\n"
          "\n"
          "  -  Fit model coefficients as an *ArrayOfGriddedField2*\n"
          "  -  Dimensions: [ n_bands ]\n"
          "\n"
          "    -   GriddedField2: [ n_band_frequencies, n_coeffs ]\n"
          "\n"
          "      -    The fit model:\n"
          "\n"
          "        -     z = p00 + p10*x + p01*y + p20*x^2\n"
          "        -     z = Xsec [m^2]\n"
          "        -     x = T / T0\n"
          "        -     y = P / P0\n"
          "        -     T0 = 1 [K]\n"
          "        -     P0 = 1 [Pa]\n"
          "        -     fitcoeffs(:, 0)           p00  [m^2]\n"
          "        -     fitcoeffs(:, 1)           p10  [m^2]\n"
          "        -     fitcoeffs(:, 2)           p01  [m^2]\n"
          "        -     fitcoeffs(:, 3)           p20  [m^2]\n"
          "\n"
          "- fitminpressures:\n"
          "\n"
          "  -  Minimum pressure available in source xsec data to generate the fit coefficients.\n"
          "  -  Dimensions: [ n_bands ]\n"
          "\n"
          "- fitmaxpressures:\n"
          "\n"
          "  -  Maximum pressure available in source xsec data to generate the fit coefficients.\n"
          "  -  Dimensions: [ n_bands ]\n"
          "\n"
          "- fitmintemperatures:\n"
          "\n"
          "  -  Minimum temperature available in source xsec data to generate the fit coefficients.\n"
          "  -  Dimensions: [ n_bands ]\n"
          "\n"
          "- fitmintemperatures:\n"
          "\n"
          "  -  Maximum temperature available in source xsec data to generate the fit coefficients.\n"
          "  -  Dimensions: [ n_bands ]\n"
          "\n"
          "fitminpressures, fitmaxpressures, fitmintemperatures and fitmaxtemperatures\n"
          "are not used to apply the model and solely serve for informational purposes.\n"),
      GROUP("ArrayOfXsecRecord")));

  wsv_data.push_back(WsvRecord(
      NAME("instrument_pol"),
      DESCRIPTION(
          "Definition of the polarisation of an instrument.\n"
          "\n"
          "The default for output is to give data for the selected Stokes\n"
          "elements (1:stokes_dim). This variable defines the polarisations\n"
          "that are actually measured, or are transmitted.\n"
          "\n"
          "The polarisation states/components are coded as\n"
          "  (0) Undefined.\n"
          "  (1) I, total intensity.\n"
          "  (2) Q, second Stokes component, Iv - Ih.\n"
          "  (3) U, third Stokes component, I+45 - I-45.\n"
          "  (4) V, forth Stokes component, Irc - Ilc\n"
          "  (5) Iv, intensity of vertically polarised component.\n"
          "  (6) Ih, intensity of horizontally polarised component.\n"
          "  (7) I+45, intensity of +45 deg linearly polarised component.\n"
          "  (8) I-45, intensity of -45 deg linearly polarised component.\n"
          "  (9) Ilhc, intensity of left-hand circularly polarised component.\n"
          "  (10) Irhc, intensity of right-hand circularly polarised component.\n"
          "\n"
          "See the documentation for definition of the Stokes vector and the\n"
          "different components.\n"
          "\n"
          "If the instrument measures, or transmits, vertical and horizontal\n"
          "components, this variable shall accordingly be set to [5,6].\n"
          "\n"
          "Conversion to Planck-BT of components 2-4 requires that component\n"
          "1 is kept, and must be included as first element.\n"
          "\n"
          "The shift from the Stokes vector can be made at any stage when of the\n"
          "sensor response set-up. The responses used must of course be adopted\n"
          "correspondingly. Or reversed, if the antenna response is defined for\n"
          "Iv or Ih it could be useful to shift polarisation as first sensor\n"
          "operation.\n"
          "\n"
          "Usage: Set by the user.\n"),
      GROUP("ArrayOfIndex")));

  wsv_data.push_back(WsvRecord(
      NAME("instrument_pol_array"),
      DESCRIPTION(
          "Multiple definition of instrument polarisation.\n"
          "\n"
          "Defined as *instrument_pol* but used when multiple polarisations\n"
          "are possible/required.\n"
          "\n"
          "Usage: Set by the user.\n"),
      GROUP("ArrayOfArrayOfIndex")));

  wsv_data.push_back(WsvRecord(
      NAME("inversion_iteration_counter"),
      DESCRIPTION(
          "Iteration counter variable for *inversion_iterate_agenda*.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("inversion_iterate_agenda"),
      DESCRIPTION(
          "Agenda recalculating spectra and Jacobian for iterative inversion methods.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("irradiance_field"),
      DESCRIPTION(
          "Irradiance field also known as flux density.\n"
          "\n"
          "Radiant flux received by a surface per unit area for each hemisphere.\n"
          "The last dimension denotes the hemispheres. The first component is\n"
          "the downward irradiance and the second component is the upward irradiance\n"
          "\n"
          "Units: W m^-2\n"
          "\n"
          "Size: [ p_grid,  lat_grid,  lon_grid,  2 ]\n"),
      GROUP("Tensor4")));

  wsv_data.push_back(
      WsvRecord(NAME("isotopologue_ratios"),
                DESCRIPTION(R"--(Contains the isotopologue ratios.

This variable is set to the default provided by *isotopologue_ratiosInitFromBuiltin*
)--"),
                GROUP("SpeciesIsotopologueRatios"),
                Species::isotopologue_ratiosInitFromBuiltin()));

  wsv_data.push_back(WsvRecord(
      NAME("iy"),
      DESCRIPTION(
          "Monochromatic pencil beam radiance spectrum.\n"
          "\n"
          "This variable holds a single spectrum, with values corresponding\n"
          "to infinite frequency and spatial resolution (compare to *y*).\n"
          "\n"
          "The variable is used to represent spectra at all positions of the\n"
          "propagation path and can e.g. temporarily hold radiation entering\n"
          "the atmosphere from space.\n"
          "\n"
          "Usage:      Used by radiative transfer methods.\n"
          "\n"
          "Unit:\n"
          " - For passive observations, as  selected by *iy_unit*.\n"
          " - For transmission calculations, same as for transmitted\n"
          "   signal.\n"
          "\n"
          "Dimensions: [ f_grid, stokes_dim ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("iyb"),
      DESCRIPTION(
          "Monochromatic pencil beam data for one measurement block.\n"
          "\n"
          "The data for all *iy* of a measurement block appended to a vector,\n"
          "following the sorting order used for *y*.\n"
          "\n"
          "Usage:      Used internally.\n"
          "\n"
          "Unit:       W / (m^2 Hz sr) or transmittance.\n"
          "\n"
          "Dimensions:\n"
          "            [ nlos * nf * stokes_dim ] where nlos is number of rows in\n"
          "            mblock_dlos, and nf is length of f_grid.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("iy_agenda_call1"),
      DESCRIPTION(
          "Flag to handle recursive calls of *iy_main_agenda*\n"
          "\n"
          "The agenda *iy_main_agenda* can be used recursively and this flag\n"
          "is used to tell the methods inside the agenda which is the primary\n"
          "call. This is handled automatically for methods using\n"
          "*iy_main_agenda*, such as *yCalc*, but the user must set this\n"
          "variable to 1 if the agenda is called directly inside the control\n"
          "file (which should be a rare case).\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("iy_aux"),
      DESCRIPTION(
          "Data auxiliary to *iy*.\n"
          "\n"
          "Different data beside the direct result of the radiative transfer\n"
          "calculations (*iy*) can be obtained by this variable. These auxilary\n"
          "data are selected by *iy_aux_vars*.\n"
          "\n"
          "Usage:      Provided by some radiative transfer methods.\n"
          "\n"
          "Dimensions: [quantity][ f_grid, stokes_dim ]\n"),
      GROUP("ArrayOfMatrix")));

  wsv_data.push_back(WsvRecord(
      NAME("iy_aux_vars"),
      DESCRIPTION(
          "Selection of quantities for *iy_aux* and when applicable also *y_aux*.\n"
          "\n"
          "Each element of this string array determines the quantity for the\n"
          "corresponding element in *iy_aux* and *y_aux* (i.e. the quantities\n"
          "are stored in the order given in *iy_aux_vars*).\n"
          "\n"
          "The possible choices vary between the methods. See the WSM you select\n"
          "for *iy_main_agenda* for the complete set of choices. Please not that\n"
          "if the calculations are done through *yCalc*, you can not select\n"
          "along-the-path variables.\n"),
      GROUP("ArrayOfString"), ArrayOfString{}));

  wsv_data.push_back(WsvRecord(
      NAME("iy_cloudbox_agenda"),
      DESCRIPTION(
          "Agenda deriving the intensity at boundary or interior of the cloudbox.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("iy_independent_beam_approx_agenda"),
      DESCRIPTION("Agenda dedicated to *iyIndependentBeamApproximation*."),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("iy_id"),
      DESCRIPTION(
          "Identification number of *iy*.\n"
          "\n"
          "This variable is intended to be an identification number for individual\n"
          "calculations of *iy*. This id-number can e.g. be used as input to \n"
          "*WriteXMLIndexed*, to link filenames to the different calculations.\n"
          "\n"
          "Some methods sets and updates *iy_id*. The general numbering scheme is::\n"
          "\n"
          "  xxxyyycba\n"
          "\n"
          "where xxx identifies the row in sensorPos/los (i.e. the mblock_index),\n"
          "yyy identifies pencil beam direction inside measurement block (should\n"
          "in general match a row in mblock_dlos), and cba identies later legs\n"
          "of total propagation paths, where a, b and c identifies secondary, tertiary\n"
          "and quaternary part, respectively. 1-based numbering is used. That is,\n"
          "the primary path of the first pencil beam of the first measurement block\n"
          "has iy_id = 001001000.\n"
          "\n"
          "Accordingly, the primary propagation path has cba = 000. If the primary path\n"
          "intersects with the surface, and the downwelling radiation is calculated\n"
          "for three directions, these secondary paths get cba = 001, 002 and 003.\n"
          "If tertiary paths appear, they have numbers such as 011. \n"
          "\n"
          "As the numbering scheme has nine positions, it is suitable to store\n"
          "files as: WriteXMLIndexed(output_file_format,iy_id,in,filename,9)\n"
          "\n"
          "Setting of *iy_id* is not yet supported together with scattering\n"
          "calculations. The value of iy_id then differs, it is either set to 0\n"
          "or keeps its value set by *yCalc*.\n"),
      GROUP("Index"), Index{0}));

  wsv_data.push_back(
      WsvRecord(NAME("iy_loop_freqs_agenda"),
                DESCRIPTION("Agenda dedicated to *iyLoopFrequencies*."),
                GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("iy_main_agenda"),
      DESCRIPTION(
          "Agenda calculating the single monochromatic pencil beam spectrum.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("iy_radar_agenda"),
      DESCRIPTION(
          "Agenda calculating pointwise backscattering.\n"),
      GROUP("Agenda")));
  
  wsv_data.push_back(WsvRecord(
      NAME("iy_space_agenda"),
      DESCRIPTION(
          "Agenda providing the downwelling radiation at the top of the atmosphere.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("iy_surface_agenda"),
      DESCRIPTION(
          "Agenda providing the upwelling radiation from the surface.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("iy_transmittance"),
      DESCRIPTION(
          "Transmittance to be included in *iy*.\n"
          "\n"
          "The calculation of *iy* can be performed over several propation path\n"
          "branches, and there can be recursive calls of *iy_main_agenda*.\n"
          "This variable gives the transmittance from the end point of the present\n"
          "branch and the sensor for such recursive cases.\n"
          "\n"
          "This variable is used purely internally. The exact usage can vary\n"
          "between different RT integration schemes.\n"
          "\n"
          "Usage:      Internally inside iy_main_agenda.\n"
          "\n"
          "Unit:       1\n"
          "\n"
          "Dimensions: [ f_grid, stokes_dim, stokes_dim ]\n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("iy_transmitter"),
      DESCRIPTION(
          "Monochromatic pencil beam radiance spectrum of transmitter signal.\n"
          "\n"
          "This variable holds a single spectrum, with values corresponding\n"
          "to infinite frequency and spatial resolution (compare to *y*).\n"
          "\n"
          "Unit:       Depend on the transmitted signal\n"
          "\n"
          "Dimensions: [ f_grid, stokes_dim ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("iy_unit"),
      DESCRIPTION(
          "Selection of output unit for radiative transfer methods.\n"
          "\n"
          "This variable allows that the unit of the output radiance/intensity\n"
          "is changed. The possible choices differ between the radiative\n"
          "methods, including not considering the variable at all.\n"
          "Accordingly, for details see the radiative method you have selected\n"
          "(e.g., *iyEmissionStandard*, *iyMC* and the like).\n"),
      GROUP("String"), String{"1"}));

  wsv_data.push_back(WsvRecord(
      NAME("iy_unit_radar"),
      DESCRIPTION(
          "Unit for radar simulations.\n"
          "\n"          
          "See the radar methods for allowed options.\n"),
      GROUP("String")));

  wsv_data.push_back(WsvRecord(
      NAME("jacobian"),
      DESCRIPTION(
          "The Jacobian matrix.\n"
          "\n"
          "The matrix holding the Jacobians of the retrieval quantities. The\n"
          "matrix has to be initialised before the retrieval quantities can be\n"
          "defined. Initialisation is done by *jacobianInit*. Retrieval quantities\n"
          "are then added with ``jacobianAdd...`` or ``retrievalAdd...`` methods.\n"
          "\n"
          "The order between rows and columns follows how data are stored in *y*\n"
          "and *x*, respectively.\n"
          "\n"
          "Units:   See the different retrieval quantities.\n"
          "\n"
          "Dimension: [ y, number of retrieval quantities and grids ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(
      WsvRecord(NAME("jacobian_agenda"),
                DESCRIPTION("Pure numerical Jacobian calculation agenda.\n"),
                GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("jacobian_do"),
      DESCRIPTION(
          "Flag to activate (clear-sky) Jacobian calculations.\n"
          "\n"
          "If this variable is set to 0, no Jacobian calculations will be done,\n"
          "even if such calculations have been set-up (through the\n"
          "jacobianAddXxx methods).\n"
          "\n"
          "Needs to be 0 if cloudy-sky (Doit) Jacobians shall be calculated.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("jacobian_quantities"),
      DESCRIPTION(
          "The retrieval quantities in the Jacobian matrix.\n"
          "\n"
          "An array of retrieval quantities for which the Jacobians are\n"
          "calculated.\n"
          "\n"
          "Usage: Quantities are added by the jacobianAdd WSMs.\n"),
      GROUP("ArrayOfRetrievalQuantity")));

  wsv_data.push_back(WsvRecord(
      NAME("jacobian_targets"),
      DESCRIPTION(
          "The partial derivatives that are computed for the Jacobian matrix.\n"
          "\n"
          "An array of jacobian targets for which the Jacobians are\n"
          "calculated.\n"
          "\n"
          "Usage: Input to absorption agendas.\n"),
      GROUP("ArrayOfJacobianTarget")));

  wsv_data.push_back(WsvRecord(NAME("lat"),
                               DESCRIPTION("A latitude.\n"
                                           "\n"
                                           "Unit:  degrees\n"),
                               GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("lat_grid"),
      DESCRIPTION(
          "The latitude grid.\n"
          "\n"
          "The latitudes for which the atmospheric fields are defined. The\n"
          "atmosphere is undefined outside the range covered by the grid.\n"
          "The grid must be sorted in increasing order, with no repetitions.\n"
          "\n"
          "Geocentric latitudes are used.\n"
          "\n"
          "For 1D calculations this vector shall be set to be empty.\n"
          "\n"
          "For 2D cases the latitudes shall be interpreted as the angular\n"
          "distance inside the orbit plane from the equator (values\n"
          "outside +-90 deg are allowed).\n"
          "\n"
          "For 3D, the valid latitude range is [-90,90].\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  degrees\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("lat_true"),
      DESCRIPTION(
          "Latitudinal geolocation for 1D and 2D data.\n"
          "\n"
          "The variables *lat_grid* and *lon_grid* contain true positions only\n"
          "for 3D. For 1D and 2D, the geographical position is given by\n"
          "*lat_true* and *lon_true*. Can be left empty when not used.\n"
          "Otherwise:\n"
          "\n"
          "   -\n"
          "     1D.\n"
          "     *lat_true* shall have length 1\n"
          "   -\n"
          "     2D.\n"
          "     Both *lat_true* and *lon_true* shall have a length matching\n"
          "     *lat_grid*. That is, *lat_true* and *lon_true* shall not be\n"
          "     seen as grids, they are vectors giving the actual lat or lon\n"
          "     for each point corresponding to *lat_grid*.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  degrees\n"),
      GROUP("Vector"), Vector{}));

  wsv_data.push_back(WsvRecord(
    NAME("lbl_checked"),
      DESCRIPTION("Flag to check if the line-by-line calculations will work\n"
                  "\n"
                  "Usage: Set manually on own risk, or use *lbl_checkedCalc*.\n"
                  "\n"
                  "Unit:  Boolean\n"),
      GROUP("Index"), Index{0}));
  
  wsv_data.push_back(
      WsvRecord(NAME("line_irradiance"),
                DESCRIPTION("Irradiance as seen by a single absorption line.\n"
                            "\n"
                            "Used internally for, e.g., NLTE effects\n"),
                GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("line_transmission"),
      DESCRIPTION("Transmission as seen by a single absorption line.\n"
                  "\n"
                  "Used internally for, e.g., NLTE effects\n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("lo"),
      DESCRIPTION(
          "The local oscillator frequency.\n"
          "\n"
          "A local oscillator frequency is used in a heterodyne system when\n"
          "the mixer folds the spectra from from radio frequencies (RF) to\n"
          "intermediate frequencies (IF).\n"
          "\n"
          "Unit:  Hz\n"
          "\n"
          "Usage: Set by the user.\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("lo_multi"),
      DESCRIPTION(
          "Local oscillator frequencies.\n"
          "\n"
          "As *lo* but describes an instrument with multiple mixers. A vector\n"
          "element for each LO. The size of this variable and\n"
          "*sideband_response_multi* shall match, and probably also\n"
          "*sideband_mode_multi*.\n"
          "\n"
          "Unit:  Hz\n"
          "\n"
          "Usage: Set by the user.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(NAME("lon"),
                               DESCRIPTION("A longitude.\n"
                                           "\n"
                                           "Unit:  degrees\n"),
                               GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("lon_grid"),
      DESCRIPTION(
          "The longitude grid.\n"
          "\n"
          "The longitudes for which the atmospheric fields are defined. The\n"
          "atmosphere is undefined outside the range covered by the grid.\n"
          "The grid must be sorted in increasing order, with no repetitions.\n"
          "\n"
          "For 1D and 2D, this WSV shall be set to be empty.\n"
          "\n"
          "Allowed values for longitudes is the range [-360,360]. The difference\n"
          "between last and first value can not exceed 360 degrees. A difference\n"
          "of exactly 360 deg. means that the complete globe is covered and no\n"
          "propagation paths will reach a longitude edge.\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  degrees\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("lon_true"),
      DESCRIPTION(
          "Longitudinal geolocation for 1D and 2D data.\n"
          "\n"
          "The variables *lat_grid* and *lon_grid* contain true positions only\n"
          "for 3D. For 1D and 2D, the geographical position is given by\n"
          "*lat_true* and *lon_true*. Can be left empty when not used.\n"
          "Otherwise:\n"
          "\n"
          "   -\n"
          "    1D.\n"
          "    *lon_true* shall have length 1\n"
          "\n"
          "   -\n"
          "     2D.\n"
          "     Both *lat_true* and *lon_true* shall have a length matching\n"
          "     *lat_grid*. That is, *lat_true* and *lon_true* shall not be\n"
          "     seen as grids, they are vectors giving the actual lat or lon\n"
          "     for each point corresponding to *lat_grid*.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  degrees\n"),
      GROUP("Vector"), Vector{}));

  wsv_data.push_back(WsvRecord(
      NAME("mag_u_field"),
      DESCRIPTION(
          "Zonal component of the magnetic field.\n"
          "\n"
          "The East-West magnetic field component. Positive values, when\n"
          "pointing eastward.\n"
          "\n"
          "Can be set to be empty, which is interpreted as zero field strength\n"
          "everywhere.\n"
          "\n"
          "Unit:       T\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ]  or [ 0 0 0 ].\n"),
      GROUP("Tensor3"), Tensor3{}));

  wsv_data.push_back(WsvRecord(
      NAME("mag_u_field_raw"),
      DESCRIPTION(
          "Raw zonal component of the magnetic field.\n"
          "\n"
          "The East-West magnetic field component. Positive values, when\n"
          "pointing eastward.\n"
          "\n"
          "Can be set to be empty, which is interpreted as zero field strength\n"
          "everywhere.\n"
          "\n"
          "Unit:       T\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ].\n"),
      GROUP("GriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("mag_v_field"),
      DESCRIPTION(
          "Meridional component of the magnetic field.\n"
          "\n"
          "The North-South magnetic field component. Positive values, when\n"
          "pointing northward.\n"
          "\n"
          "Can be set to be empty, which is interpreted as zero field strength\n"
          "everywhere.\n"
          "\n"
          "Unit:       T\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ]  or [ 0 0 0 ].\n"),
      GROUP("Tensor3"), Tensor3{}));

  wsv_data.push_back(WsvRecord(
      NAME("mag_v_field_raw"),
      DESCRIPTION(
          "Raw meridional component of the magnetic field.\n"
          "\n"
          "The North-South magnetic field component. Positive values, when\n"
          "pointing northward.\n"
          "\n"
          "Can be set to be empty, which is interpreted as zero field strength\n"
          "everywhere.\n"
          "\n"
          "Unit:       T\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ].\n"),
      GROUP("GriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("mag_w_field"),
      DESCRIPTION(
          "Vertical component of the magnetic field.\n"
          "\n"
          "Positive values, when pointing upward.\n"
          "\n"
          "Can be set to be empty, which is interpreted as zero field strength\n"
          "everywhere.\n"
          "\n"
          "Unit:       T\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ]  or [ 0 0 0 ].\n"),
      GROUP("Tensor3"), Tensor3{}));

  wsv_data.push_back(WsvRecord(
      NAME("mag_w_field_raw"),
      DESCRIPTION(
          "Raw vertical component of the magnetic field.\n"
          "\n"
          "Positive values, when pointing upward.\n"
          "\n"
          "Can be set to be empty, which is interpreted as zero field strength\n"
          "everywhere.\n"
          "\n"
          "Unit:       T\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ].\n"),
      GROUP("GriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("main_agenda"),
      DESCRIPTION("Agenda corresponding to the entire controlfile.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("mblock_dlos"),
      DESCRIPTION(
          "The set of angular pencil beam directions for each measurement block.\n"
          "\n"
          "The relative angles in this variable are angular off-sets with\n"
          "respect to the angles in *sensor_los*. Defined as *dlos* but is\n"
          "allowed to only have a single column, as described below.\n"
          "\n"
          "The first column holds the relative zenith angle. This column is\n"
          "mandatory for all atmospheric dimensionalities. For 3D, there can\n"
          "also be a second column, giving relative azimuth angles. If this\n"
          "column is not present (for 3D) zero azimuth off-sets are assumed.\n"
          "\n"
          "Usage: Set by the user or output of antenna WSMs.\n"
          "\n"
          "Unit:  degrees\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("mblock_index"),
      DESCRIPTION(
          "Measurement block index. \n"
          "\n"
          "Used to tell agendas the index of present measurement block.\n"
          "\n"
          "Usage: Used internally.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_antenna"),
      DESCRIPTION(
          "Antenna pattern description for dedicated MC calculaions.\n"
          "\n"
          "Usage: Input to MCGeneral. Set by *mc_antennaSetGaussian* and similar methods.\n"),
      GROUP("MCAntenna")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_error"),
      DESCRIPTION("Error in simulated *y* when using a Monte Carlo approach.\n"
                  "\n"
                  "Usage: Output from Monte Carlo functions. \n"
                  "\n"
                  "Units: Depends on *iy_unit*.\n"
                  "\n"
                  "Size:  [ stokes_dim ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_iteration_count"),
      DESCRIPTION(
          "Counts the number of iterations (or photons) used in the MC\n"
          "scattering algorithm.\n"
          "\n"
          "Usage: Set by MCGeneral and other MC methods.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_max_iter"),
      DESCRIPTION("The maximum number of iterations allowed for Monte Carlo\n"
                  "calculations.\n"
                  "\n"
                  "Usage: Set by the user.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_max_scatorder"),
      DESCRIPTION("The maximum scattering order allowed for Monte Carlo\n"
                  "radar calculations.\n"
                  "\n"
                  "Usage: Set by the user.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_max_time"),
      DESCRIPTION("The maximum time allowed for Monte Carlo calculations.\n"
                  "\n"
                  "Usage: Set by the user.\n"
                  "\n"
                  "Unit: s\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_min_iter"),
      DESCRIPTION("The minimum number of iterations allowed for Monte Carlo\n"
                  "calculations.\n"
                  "\n"
                  "Usage: Set by the user.\n"),
      GROUP("Index"), Index{100}));

  wsv_data.push_back(
      WsvRecord(NAME("mc_points"),
                DESCRIPTION(
                    //FIXMEDOC
                    "Source to emission, position.\n"
                    "\n"
                    "Counts the number of MC endpoints in each grid cell.\n"
                    "\n"
                    "Usage: Set by MCGeneral and other MC methods.\n"),
                GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_scat_order"),
      DESCRIPTION(
          "Number of atmospheric scattering events between emission point and sensor.\n"
          "\n"
          "The first element gives the number of cases with zero scattering events,\n"
          "the second the number of single scattering cases etc.\n"
          "\n"
          "Scattering orders above what the variable can hold are not stored at all.\n"
          "The number of such cases can be determined by comparing\n"
          "*mc_iteration_count* with the sum of the elements in this array.\n"
          "\n"
          "Usage: Set by MCGeneral and other MC methods.\n"),
      GROUP("ArrayOfIndex")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_source_domain"),
      DESCRIPTION(
          "Rough classification of source to emission.\n"
          "\n"
          "This is an array of length 4, where the elements in order represent\n"
          "space, the surface, atmospheric gas and atmospheric particle.\n"
          "The distinction between the two last elements is if the emission\n"
          "is associated with *vmr_field* or *pnd_field*.\n"
          "\n"
          "The values of the array give the number of cases where the emission\n"
          "source was found to be inside each \"domain\".\n"
          "\n"
          "Usage: Set by MCGeneral and other MC methods.\n"),
      GROUP("ArrayOfIndex")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_seed"),
      DESCRIPTION("The integer seed for the random number generator used by\n"
                  "Monte Carlo methods.\n"
                  "\n"
                  "Usage: Set by MCSetSeed.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_std_err"),
      DESCRIPTION(
          "Target precision (1 std. dev.) for Monte Carlo calculations.\n"
          "\n"
          "Usage: Set by the user.\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_y_tx"),
      DESCRIPTION("Normalized Stokes vector for transmittance (e.g., radar).\n"
                  "\n"
                  "The first element (intensity) should have a value of 1."
                  "\n"
                  "Usage: Set by user. \n"
                  "\n"
                  "Units: Unitless.\n"
                  "\n"
                  "Size:  [ stokes_dim ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("mc_taustep_limit"),
      DESCRIPTION(
          "Defines an upper step length in terms of optical thickness for Monte "
          "Carlo calculations.\n"
          "\n"
          "Usage: Set by the user.\n"),
      GROUP("Numeric"), Numeric{0.1}));

  wsv_data.push_back(WsvRecord(
      NAME("met_amsu_data"),
      DESCRIPTION(
          "The AMSU data set.\n"
          "\n"
          "This is intended as input for the method ybatchMetProfiles. It holds the\n"
          "latitude, longitude, satellite zenith angle and amsu-b corrected and \n"
          "uncorrected brightness temperatures.  It also has information about \n"
          "the particular pixel corresponds to a land or sea point.  This will be \n"
          "read in the method ybatchMetProfiles and the profiles corresponding to \n"
          "each latitude and longitude will be read in.\n"
          "\n"
          "See documentation of WSM *ybatchMetProfiles* for more information.\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("met_mm_antenna"),
      DESCRIPTION(
          "The antenna beam width for meteorological millimeter instruments.\n"
          "\n"
          "This Vector must match the number and order of channels in\n"
          "*met_mm_backend*.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  [ Hz ]\n"
          "\n"
          "Size:  [ number of channels ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("met_mm_backend"),
      DESCRIPTION(
          "Backend description for meteorological millimeter sensors with passbands.\n"
          "\n"
          "This is a compact description of a passband-type sensor, e.g. AMSU-A. The matrix\n"
          "contains one row for each instrument channel. Each row contains four elements:\n\n"
          " - LO position [Hz]\n"
          " - first offset from the LO [Hz]\n"
          " - second offset from the LO+offset1 [Hz]\n"
          " - channel width [Hz]\n"
          "\n"
          "Overview::\n\n"
          "                            LO\n"
          "                             |\n"
          "                 offset1     |    offset1\n"
          "             ----------------+----------------\n"
          "             |                               |\n"
          "             |                               |\n"
          "    offset2  |  offset2             offset2  |  offset2\n"
          "    ---------+---------             ---------+---------\n"
          "    |                 |             |                 |\n"
          "    |                 |             |                 |\n"
          "  #####             #####         #####             #####\n"
          "  width             width         width             width\n"
          "\n"
          "For a sensor with 1 passband, offset1 and offset2 are zero.\n"
          "For a sensor with 2 passbands, only offset2 is zero.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit: All entries in Hz.\n"
          "\n"
          "Size: [number of channels, 4]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("met_mm_polarisation"),
      DESCRIPTION(
          "The polarisation for meteorological millimeter instruments.\n"
          "\n"
          "This array must match the number and order of channels in\n"
          "*met_mm_backend*.\n"
          "\n"
          "Possible values:\n\n"
          "- ``\"V\"``: Vertical polarisation\n"
          "- ``\"H\"``: Horizontal polarisation\n"
          "- ``\"LHC\"``: Left-hand circular polarisation\n"
          "- ``\"RHC\"``: Right-hand circular polarisation\n"
          "- ``\"AMSU-V\"``: Vertical polarisation dependening on AMSU zenith angle\n"
          "- ``\"AMSU-H\"``: Horizontal polarisation dependening on AMSU zenith angle\n"
          "- ``\"ISMAR-V\"``: Vertical polarisation dependening on ISMAR zenith angle\n"
          "- ``\"ISMAR-H\"``: Horizontal polarisation dependening on AMSU zenith angle\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  [ String ]\n"
          "\n"
          "Size:  [ number of channels ]\n"),
      GROUP("ArrayOfString")));

  wsv_data.push_back(
      WsvRecord(NAME("met_profile_calc_agenda"),
                DESCRIPTION("Agenda for metoffice profile calculations.\n"),
                GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("level0_data"),
      DESCRIPTION("List of L0 data.  Can be of any type.\n"
                  "It is method-dependent how this is used to calibrated to L1\n"),
      GROUP("ArrayOfVector")));

  wsv_data.push_back(WsvRecord(
      NAME("level0_time"),
      DESCRIPTION("List of L0 times.  Should be in UTC.\n"
                  "It is method-dependent how this is used to calibrated to L1\n"),
      GROUP("ArrayOfTime")));

  wsv_data.push_back(WsvRecord(
      NAME("lm_ga_history"),
      DESCRIPTION(
          "The series of gamma values for a Marquardt-levenberg inversion.\n"
          "\n"
          "The values are stored following iteration order, i.e. the first\n"
          "is the gamma factor for the first iteration etc.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("molarmass_dry_air"),
      DESCRIPTION(
          "The average molar mass of dry air.\n"
          "\n"
          "This could also be referred to as the average molecular weight for\n"
          "dry air. The definition of \"dry air\" can differ between planets and\n"
          "methods using the WSV. For Earth, this should be a value around\n"
          "28.97.\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("nlte_level_identifiers"),
      DESCRIPTION("An array of non-lte quantum identifiers for levels matching\n"
                  "*nlte_field_raw* and on request *nlte_vibrational_energies*.\n"),
      GROUP("ArrayOfQuantumIdentifier"), ArrayOfQuantumIdentifier{}));

  wsv_data.push_back(WsvRecord(
      NAME("nlte_vibrational_energies"),
      DESCRIPTION("An list of vibrational energies matching\n"
                  "*nlte_level_identifiers* and *nlte_field_raw* or being 0.\n"),
      GROUP("Vector"))); 

  wsv_data.push_back(WsvRecord(
      NAME("collision_line_identifiers"),
      DESCRIPTION(
          //FIXMEDOC
          "An array of quantum identifiers for finding collisional rates\n"
          "in *collision_coefficients*\n"),
      GROUP("ArrayOfQuantumIdentifier")));

  wsv_data.push_back(
      WsvRecord(NAME("collision_coefficients"),
                DESCRIPTION(
                    //FIXMEDOC
                    "An array of coefficients for effective collisions\n"),
                GROUP("ArrayOfArrayOfGriddedField1")));

  wsv_data.push_back(
      WsvRecord(NAME("nelem"),
                DESCRIPTION("Number of elements of a Vector or Array.\n"),
                GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("ncols"),
      DESCRIPTION(
          "Number of columns (elements in lowest dimension) of a Matrix or Tensor.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("nrows"),
      DESCRIPTION(
          "Number of rows (elements in 2nd lowest dimension) of a Matrix or Tensor.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("npages"),
      DESCRIPTION("Number of elements in 3rd lowest dimension of a Tensor.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("nbooks"),
      DESCRIPTION("Number of elements in 4th lowest dimension of a Tensor.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("nshelves"),
      DESCRIPTION("Number of elements in 5th lowest dimension of a Tensor.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("nvitrines"),
      DESCRIPTION("Number of elements in 6th lowest dimension of a Tensor.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("nlibraries"),
      DESCRIPTION("Number of elements in 7th lowest dimension of a Tensor.\n"),
      GROUP("Index")));

  wsv_data.push_back(
      WsvRecord(NAME("nlte_do"),
                DESCRIPTION("Flag to perform Non-LTE calculations.\n"),
                GROUP("Index"), Index{0}));

  wsv_data.push_back(WsvRecord(
      NAME("nlte_source"),
      DESCRIPTION(
          "Variable to contain the additional source function due to NLTE effects.\n"
          "\n"
          "Dimensions: [nza, naa, nf, stokes_dim]\n"),
      GROUP("StokesVector")));

  wsv_data.push_back(WsvRecord(
      NAME("oem_diagnostics"),
      DESCRIPTION(
          "Basic diagnostics of an OEM type inversion.\n"
          "\n"
          "This is a vector of length 5, having the elements (0-based index):\n"
          "  0. Convergence status, with coding\n"
          "       - 0 = converged\n"
          "       - 1 = max iterations reached\n"
          "       - 2 = max gamma of LM reached\n"
          "       - 9 = some error when calling *inversion_iterate_agenda*\n"
          "       - 99 = too high start cost.\n"
          "  1. Start value of cost function.\n"
          "  2. End value of cost function.\n"
          "  3. End value of y-part of cost function.\n"
          "  4. Number of iterations used.\n"
          "\n"
          "See WSM *OEM* for a definition of \"cost\". Values not calculated\n"
          "are set to NaN.\n"),
      GROUP("Vector")));
  wsv_data.push_back(
      WsvRecord(NAME("oem_errors"),
                DESCRIPTION("Errors encountered during OEM execution.\n"),
                GROUP("ArrayOfString")));

  wsv_data.push_back(WsvRecord(
      NAME("output_file_format"),
      DESCRIPTION(
          "Output file format.\n"
          "\n"
          "This variable sets the format for output files. It could be set to:\n\n"
          "- \"ascii\": for plain xml files\n"
          "- \"zascii\": for zipped xml files\n"
          "- \"binary\": for binary.\n"
          "\n"
          "To change the value of this variable use the workspace methods\n"
          "*output_file_formatSetAscii*, *output_file_formatSetZippedAscii*, and\n"
          "*output_file_formatSetBinary*\n"),
      GROUP("String"), String{"ascii"}));

  wsv_data.push_back(WsvRecord(
      NAME("particle_bulkprop_field"),
      DESCRIPTION(
          "Container for various data that describes scattering bulk properties.\n"
          "\n"
          "The number and order of bulk properties is free, as long as the data are\n"
          "consistent with the content of *particle_bulkprop_names*. \n"
          "\n"
          "The data shall be given on the standard atmospheric grids. When actually\n"
          "used, this variable must have zeros at all positions outside and at the\n"
          "border of the cloudbox.\n"
          "\n"
          "Dimensions: [ particle_bulkprop_names, p_grid, lat_grid, lon_grid ]\n"),
      GROUP("Tensor4"), Tensor4{}));

  wsv_data.push_back(WsvRecord(
      NAME("particle_bulkprop_names"),
      DESCRIPTION(
          "Identification of the data in *particle_bulkprop_field*.\n"
          "\n"
          "This variable assigns a name to each field in *particle_bulkprop_field*.\n"
          "The naming is totally free. If two fields are given the same name, the\n"
          "first one will be selected.\n"
          "\n"
          "Dimensions: length should match book-dimension of *particle_bulkprop_field*\n"),
      GROUP("ArrayOfString"), ArrayOfString{}));

  wsv_data.push_back(WsvRecord(
      NAME("particle_masses"),
      DESCRIPTION(
          "The mass of individual particles (or bulks).\n"
          "\n"
          "Each row corresponds to a scattering element (i.e. an element in\n"
          "*scat_data*). The user is free to define different mass\n"
          "categories and assign a mass for each category. Each column\n"
          "of *particle_masses* corresponds to such a mass category. A scattering\n"
          "element can have a non-zero mass for more than one category.\n"
          "\n"
          "For example, if you work with clouds, your mass categories could\n"
          "be ice and liquid, corresponding to IWC and LWC, respectively.\n"
          "The mass of particles inside the melting layer, having a mixed\n"
          "phase, could be divided between the two columns of the matrix.\n"
          "\n"
          "Shall either be empty, or have a row size consistent with the\n"
          "scattering variables (*scat_data*, *pnd_field*).\n"
          "\n"
          "Usage:      Set by the user.\n"
          "\n"
          "Unit:       kg\n"
          "\n"
          "Dimensions: [number of scattering elements, number of mass categories]\n"),
      GROUP("Matrix"), Matrix{}));

  wsv_data.push_back(WsvRecord(
      NAME("pha_mat"),
      DESCRIPTION(
          "Ensemble averaged phase matrix.\n"
          "\n"
          "This workspace variable represents the actual physical phase\n"
          "matrix (averaged over all scattering elements) for given propagation\n"
          "directions. It is calculated in the method *pha_matCalc*.\n"
          "\n"
          "See ARTS user guide (AUG) for further information. Use the index to find\n"
          "where this variable is discussed. The variable is listed as a subentry\n"
          "to \"workspace variables\".\n"
          "\n"
          "Usage:      Output of the method *pha_matCalc*\n"
          "\n"
          "Unit:        m^2\n"  //FIXME: really m2? not 1/m?
          "\n"
          "Dimensions: [za_grid, aa_grid, stokes_dim, stokes_dim]\n"),
      GROUP("Tensor4")));

  wsv_data.push_back(WsvRecord(
      NAME("pha_mat_doit"),
      DESCRIPTION(
          "Ensemble averaged phase matrix for DOIT calculation.\n"
          "\n"
          "This workspace variable represents the actual physical phase\n"
          "matrix (averaged over all scattering elements) for given incident and \n"
          "propagation directions. It is calculated in the method *DoitScatteringDataPrepare*.\n"
          "\n"
          "See ARTS user guide (AUG) for further information."
          "\n"
          "Usage:      Output of the method *pha_matCalc*\n"
          "\n"
          "Unit:        m^2\n"  //FIXME: really m2? not 1/m?
          "\n"
          "Dimensions: [T,za_grid, aa_grid, za_grid, aa_grid, stokes_dim, stokes_dim]\n"),
      GROUP("Tensor7")));

  wsv_data.push_back(WsvRecord(
      NAME("pha_mat_spt"),
      DESCRIPTION(
          "Phase matrix for all individual scattering elements.\n"
          "\n"
          "This variable contains the elements of phase matrix for all individual\n"
          "scattering elements for given propagation directions. It is the\n"
          "calculated in the agenda *pha_mat_spt_agenda*. The elements of the phase\n"
          "matrix are calculated from the single scattering data.\n"
          "\n"
          "See ARTS user guide (AUG) for further information.\n"
          "\n"
          "Usage:      Input and Output of the pha_mat_sptFrom* methods\n"
          "\n"
          "Unit:       m^2\n"  //FIXME: really m2? not 1/m?
          "\n"
          "Dimensions: [number of scattering elements, za_grid, aa_grid, stokes_dim, stokes_dim]\n"),
      GROUP("Tensor5")));

  wsv_data.push_back(WsvRecord(
      NAME("pha_mat_spt_agenda"),
      DESCRIPTION(
          "Agenda calculates the phase matrix for individual scattering elements.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("pha_mat_sptDOITOpt"),
      DESCRIPTION(
          "Interpolated phase matrix.\n"
          "\n"
          "This variable contains the data of the phase matrix in the \n"
          "scattering frame interpolated on the actual frequency (the variable\n"
          "is used inside *doit_mono_agenda*) and also interpolated on all \n"
          "possible scattering angles following from all combinations of \n"
          "*za_grid* and *aa_grid*. \n"
          "\n"
          "Usage:      Input of the method *pha_mat_sptFromDataDOITOpt*\n"
          "\n"
          "Unit:        m^2\n"  //FIXME: really m2? not 1/m?
          "\n"
          "Dimensions: \n"
          "\n"
          "- [number of scattering elements]\n"
          "\n"
          "  - [T, za_grid, aa_grid, za_grid, aa_grid, stokes_dim, stokes_dim]\n"),
      GROUP("ArrayOfTensor7")));

  wsv_data.push_back(WsvRecord(
      NAME("planet_rotation_period"),
      DESCRIPTION(
          "The sidereal rotation period of the planet.\n"
          "\n"
          "This is time that it takes for the planet to complete one revolution\n"
          "around its axis of rotation relative to the suns. For Earth, this\n"
          "is a value roughly 4 min less than 24 h.\n"
          "\n"
          "A negative value signifies a retrograde rotation, i.e. opposite to\n"
          "the rotation of Earth.\n\n"
          "Unit:   s\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("pnd_agenda_array"),
      DESCRIPTION(
          "Mapping of particle bulk properties to number density data.\n"
          "\n"
          "The length of this agenda array shall match the size of *scat_species*.\n"
          "That is there is a \"pnd-agenda\" associated with each scattering species.\n"
          "\n"
          "In short, each agenda takes some bulk property data as input, and returns\n"
          "particle number densities for all scattering elements of the species.\n"
          "See further *pnd_agenda_input* and associated variables.\n"),
      GROUP("ArrayOfAgenda")));

  wsv_data.push_back(WsvRecord(
      NAME("pnd_agenda_input"),
      DESCRIPTION(
          "The variable input to one element of *pnd_agenda_array*.\n"
          "\n"
          "The column dimension corresponds to the input to the underlying\n"
          "particle size distribution method. For example, the first column\n"
          "can hold ice water content values, and the second one temperature\n"
          "data.\n"
          "\n"
          "Temperatures are handled by *pnd_agenda_input_t* and shall not be\n"
          "included in this variable.\n"
          "\n"
          "Each row corresponds to a position. That is, the methods in the\n"
          "pnd-agendas are expected to process multiple points in one call.\n"
          "\n"
          "Dimensions: [ n_points, n_input_variables ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("pnd_agenda_input_t"),
      DESCRIPTION(
          "Temperature input to one element of *pnd_agenda_array*.\n"
          "\n"
          "This WSV works as *pnd_agenda_input* but holds a specific quantity,\n"
          "temperature.\n"
          "\n"
          "Each element corresponds to a position. That is, the methods in the\n"
          "pnd-agendas are expected to process multiple points in one call.\n"
          "\n"
          "Dimensions: [ n_points ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("pnd_agenda_array_input_names"),
      DESCRIPTION(
          "Naming of all input expected by *pnd_agenda_array*.\n"
          "\n"
          "This variable contains *pnd_agenda_input_names* for each agenda\n"
          "element in *pnd_agenda_array*.\n"
          "\n"
          "Dimension: [ n_scattering_species ][ n_input_variables ]\n"),
      GROUP("ArrayOfArrayOfString")));

  wsv_data.push_back(WsvRecord(
      NAME("pnd_agenda_input_names"),
      DESCRIPTION(
          "Naming of (existing or expected) data in *pnd_agenda_input*.\n"
          "\n"
          "The strings of this variable refer to the corresponding column in\n"
          "*pnd_agenda_input*.\n"
          "\n"
          "Dimension: [ n_input_variables ]\n"),
      GROUP("ArrayOfString")));

  wsv_data.push_back(WsvRecord(
      NAME("pnd_data"),
      DESCRIPTION(
          "Particle number density values for a set of points.\n"
          "\n"
          "The variable contains particle number density data for one scattering\n"
          "species. The row dimension corresponds to different positions, in the\n"
          "same way as *pnd_agenda_input* is defined.\n"
          "\n"
          "Dimensions: [ n_points, n_scattering_elements ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("pnd_field"),
      DESCRIPTION(
          "Particle number density field.\n"
          "\n"
          "This variable holds the particle number density fields for all\n"
          "scattering elements being read in the WSMs\n"
          "*ScatElementsPndAndScatAdd* or *ScatSpeciesPndAndScatAdd* and\n"
          "interpolated to the calculation grids *p_grid*, *lat_grid*, and\n"
          "*lon_grid* inside the cloudbox. An alternative method to create\n"
          "*pnd_field* is *pnd_fieldCalcFromParticleBulkProps*.\n"
          "\n"
          "Total number and order of scattering elements in *pnd_field* and (the\n"
          "flattened) *scat_data* has to be identical.\n"
          "\n"
          "Note: To ensure that no particles exist outside the cloudbox,\n"
          "*pnd_field* is required to be 0 at its outer limits (corresponding\n"
          "to the *cloudbox_limits*).\n"
          "\n"
          "Usage:      Set by user or output of *pnd_fieldCalcFromParticleBulkProps*\n"
          "\n"
          "Unit:        m^-3\n"
          "\n"
          "Size:\n"
          " [number of scattering elements, \n"
          " ( *cloudbox_limits* [1] - *cloudbox_limits* [0]) +1, \n"
          " ( *cloudbox_limits* [3] - *cloudbox_limits* [2]) +1, \n"
          " ( *cloudbox_limits* [5] - *cloudbox_limits* [4]) +1 ] \n"),
      GROUP("Tensor4")));

  wsv_data.push_back(WsvRecord(
      NAME("pnd_size_grid"),
      DESCRIPTION(
          "The particle sizes associated with *pnd_data*.\n"
          "\n"
          "This variable holds the size of each scattering element considered.\n"
          "Size can be defined differently, depending on particle size distribution\n"
          "used. Most common choices should by equivalent diameter, maximum diameter\n"
          "and mass.\n"
          "\n"
          "Dimension: [ n_sizes ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("pnd_field_raw"),
      DESCRIPTION(
          "The particle number density field raw data.\n"
          "\n"
          "This variable contains the particle number density data for all\n"
          "considered scattering elements. *pnd_field_raw* is an Array of\n"
          "GriddedField3. It includes one GriddedField3 for each scattering\n"
          "element, which contains both the data and the corresponding grids.\n"
          "\n"
          "Usage:\n"
          "  Set by the user. Input to methods *ScatElementsPndAndScatAdd* and \n"
          "  *ScatSpeciesPndAndScatAdd*\n"
          "\n"
          "Unit:  m^-3\n"
          "\n"
          "Size:\n"
          "\n"
          "- Array[number of scattering elementst]\n"
          "\n"
          "  - GriddedField3\n"
          "  - [number of pressure levels]\n"
          "  - [number of latitudes]\n"
          "  - [number of longitudes]\n"
          "  - [number of pressure levels, number of latitudes, number of longitudes]\n"),
      GROUP("ArrayOfGriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("ppath"),
      DESCRIPTION(
          "The propagation path for one line-of-sight.\n"
          "\n"
          "This variable describes the total (pencil beam) propagation path for\n"
          "a given combination of starting point and line-of-sight. The path is\n"
          "described by a data structure of type Ppath. This structure contains\n"
          "also additional fields to faciliate the calculation of spectra and\n"
          "interpolation of the atmospheric fields.\n"
          "\n"
          "The data struture is too extensive to be described here, but it is\n"
          "described carefully in the ARTS user guide (AUG). Use the index to\n"
          "find where the data structure, Ppath, for propagation paths is \n"
          "discussed. It is listed as a subentry to \"data structures\".\n"
          "\n"
          "Usage: Output from *ppath_agenda*.\n"),
      GROUP("Ppath")));

  wsv_data.push_back(
      WsvRecord(NAME("ppath_agenda"),
                DESCRIPTION("Agenda calculating complete propagation paths.\n"),
                GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("ppath_field"),
      DESCRIPTION(
          "An array meant to build up the necessary geometries for radiative\n"
          "field calculations.\n"
          "\n"
          "Can be ordered or not\n"
          "\n"
          "Size: user-defined\n"),
      GROUP("ArrayOfPpath")));

  wsv_data.push_back(WsvRecord(
      NAME("ppath_inside_cloudbox_do"),
      DESCRIPTION(
          "Flag to perform ray tracing inside the cloudbox.\n"
          "\n"
          "Standard propagation path calculations stop at the boundary of the\n"
          "cloudbox, or stop directly if started inside the cloudbox. This WSV\n"
          "allows scattering methods to obtain propagation paths inside the\n"
          "cloudbox. Hence, this variable is for internal usage primarily.\n"
          "\n"
          "Usage: For communication between modules of arts.\n"),
      GROUP("Index"), Index{0}));

  wsv_data.push_back(WsvRecord(
      NAME("ppath_lmax"),
      DESCRIPTION(
          "Maximum length between points describing propagation paths.\n"
          "\n"
          "See *ppath_stepGeometric* for a description of this variable.\n"
          "\n"
          "Usage: Ppath methods such as *ppath_stepGeometric*.\n"),
      GROUP("Numeric"), Numeric{10e3}));

  wsv_data.push_back(WsvRecord(
      NAME("ppath_lraytrace"),
      DESCRIPTION(
          "Maximum length of ray tracing steps when determining propagation\n"
          "paths.\n"
          "\n"
          "See *ppath_stepRefractionBasic* for a description of this variable.\n"
          "\n"
          "Usage: Refraction ppath methods such as *ppath_stepRefractionBasic*.\n"),
      GROUP("Numeric"), Numeric{1e3}));

  wsv_data.push_back(WsvRecord(
      NAME("ppath_step"),
      DESCRIPTION(
          "A propagation path step.\n"
          "\n"
          "The main intention of this variable is communication with the agenda\n"
          "*ppath_step_agenda*.\n"
          "\n"
          "See *ppath_step_agenda* for more information on this variable and\n"
          "the calculation of propagation paths. Or read the chapter on\n"
          "propagation paths in the ARTS user guide.\n"
          "\n"
          "Usage:   In/output to/from *ppath_step_agenda*.\n"
          "\n"
          "Members: See AUG.\n"),
      GROUP("Ppath")));

  wsv_data.push_back(
      WsvRecord(NAME("ppath_step_agenda"),
                DESCRIPTION("Agenda calculating a propagation path step.\n"),
                GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("ppvar_f"),
      DESCRIPTION(
          "Doppler adjusted frequencies along the propagation path.\n"
          "\n"
          "See *ppvar_p* for a general description of WSVs of ppvar-type.\n"
          "\n"
          "Dimension: [ number of frequencies, ppath.np ]\n"
          "\n"
          "Usage: Output of radiative transfer methods.\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("ppvar_iy"),
      DESCRIPTION(
          "iy-values along the propagation path.\n"
          "\n"
          "See *ppvar_p* for a general description of WSVs of ppvar-type.\n"
          "\n"
          "Dimension: [ number of frequencies, stokes_dim, ppath.np ]\n"
          "\n"
          "Usage: Output of radiative transfer methods.\n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("ppvar_mag"),
      DESCRIPTION(
          "Magnetic field along the propagation path.\n"
          "\n"
          "See *ppvar_p* for a general description of WSVs of ppvar-type.\n"
          "\n"
          "Dimension: [ 3, ppath.np ]\n"
          "\n"
          "Usage: Output of radiative transfer methods.\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("ppvar_nlte"),
      DESCRIPTION(
          "Non-LTE temperatures/ratios along the propagation path.\n"
          "\n"
          "See *ppvar_p* for a general description of WSVs of ppvar-type.\n"
          "\n"
          "Dimension: [ number of non-lte temperatures, 1, 1, ppath.np ]\n"
          "\n"
          "Usage: Output of radiative transfer methods.\n"),
      GROUP("EnergyLevelMap")));

  wsv_data.push_back(WsvRecord(
      NAME("ppvar_p"),
      DESCRIPTION(
          "Pressure along the propagation path.\n"
          "\n"
          "ppvar stands for propagation path variable. The variables named in is\n"
          "way describe the atmosphere and its properties at each point of the\n"
          "propagation path\n"
          "\n"
          "Dimension: [ ppath.np ]\n"
          "\n"
          "Usage: Output of radiative transfer methods.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("ppvar_pnd"),
      DESCRIPTION(
          "PND values along the propagation path.\n"
          "\n"
          "See *ppvar_p* for a general description of WSVs of ppvar-type.\n"
          "\n"
          "Dimension: [ number of scattering elements, ppath.np ]\n"
          "\n"
          "Usage: Output of radiative transfer methods.\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("ppvar_optical_depth"),
      DESCRIPTION(
          "The optical depth between the sensor and each point of the propagation path.\n"
          "\n"
          "Returned as the one-way optical depth even in the case of radar\n"
          "simulations. Just a scalar value, i.e. no polarisation information is\n"
          "provided.\n"
          "\n"
          "See *ppvar_p* for a general description of WSVs of ppvar-type.\n"
          "\n"
          "Dimension: [ ppath.np, f_grid]\n"
          "\n"
          "Usage: Output of radiative transfer methods.\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("ppvar_t"),
      DESCRIPTION(
          "Temperature along the propagation path.\n"
          "\n"
          "See *ppvar_p* for a general description of WSVs of ppvar-type.\n"
          "\n"
          "Dimension: [ ppath.np ]\n"
          "\n"
          "Usage: Output of radiative transfer methods.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("ppvar_trans_cumulat"),
      DESCRIPTION(
          "The transmittance between the sensor and each point of the propagation path.\n"
          "\n"
          "See *ppvar_p* for a general description of WSVs of ppvar-type.\n"
          "\n"
          "Dimension: [ ppath.np, f_grid, stokes_dim, stokes_dim ]\n"
          "\n"
          "Usage: Output of radiative transfer methods.\n"),
      GROUP("Tensor4")));

  wsv_data.push_back(WsvRecord(
      NAME("ppvar_trans_partial"),
      DESCRIPTION(
          "The transmittance between the points along the propagation path.\n"
          "\n"
          "To maintain consistency in size also this variable stores np transmissivities,\n"
          "while there are only np-1 distances between the points of the ppath. The\n"
          "extra values placed at index 0 and can be seen as the transmissivities\n"
          "between the sensor and the start of the ppath. These transmissivities\n"
          "are always unity. That is, the transmissivities between ppath point i and i+1\n"
          "are found at index i+1 in *ppvar_trans_partial*.\n"
          "\n"
          "See *ppvar_p* for a general description of WSVs of ppvar-type.\n"
          "\n"
          "Dimension: [ ppath.np, f_grid, stokes_dim, stokes_dim ]\n"
          "\n"
          "Usage: Output of radiative transfer methods.\n"),
      GROUP("Tensor4")));

  wsv_data.push_back(WsvRecord(
      NAME("ppvar_vmr"),
      DESCRIPTION(
          "VMR values along the propagation path.\n"
          "\n"
          "See *ppvar_p* for a general description of WSVs of ppvar-type.\n"
          "\n"
          "Dimension: [ number of abs. species, ppath.np ]\n"
          "\n"
          "Usage: Output of radiative transfer methods.\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("ppvar_wind"),
      DESCRIPTION(
          "Winds along the propagation path.\n"
          "\n"
          "See *ppvar_p* for a general description of WSVs of ppvar-type.\n"
          "\n"
          "Dimension: [ 3, ppath.np ]\n"
          "\n"
          "Usage: Output of radiative transfer methods.\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("predefined_model_data"),
      DESCRIPTION(
          R"--(This contains predefined model data.

Can currently only contain data for new MT CKD models of water.
)--"),
      GROUP("PredefinedModelData"), PredefinedModelData{}));

  wsv_data.push_back(WsvRecord(
      NAME("spectral_radiance_profile_operator"),
      DESCRIPTION(
          R"--(An operator to create a spectral radiance profile.

This is an experimental solution.
)--"),
      GROUP("SpectralRadianceProfileOperator")));

  wsv_data.push_back(WsvRecord(
      NAME("propmat_clearsky"),
      DESCRIPTION(
          "This contains the absorption coefficients for one point in the atmosphere.\n"
          "\n"
          "Dimensions: [naa, nza, nf, f(stokes_dim)]\n"
          "\n"
          "Unit: 1/m\n"),
      GROUP("PropagationMatrix")));

  wsv_data.push_back(
      WsvRecord(NAME("propmat_clearsky_agenda_checked"),
                DESCRIPTION("OK-flag for *propmat_clearsky_agenda*.\n"
                            "\n"
                            "Set by *propmat_clearsky_agenda_checkedCalc*.\n"),
                GROUP("Index"), Index{0}));

  wsv_data.push_back(WsvRecord(
      NAME("propmat_clearsky_agenda"),
      DESCRIPTION("Agenda calculating the absorption coefficient matrices.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("propmat_clearsky_field"),
      DESCRIPTION(
          "Gas absorption field.\n"
          "\n"
          "Contains the (polarized) gas absorption coefficients for all species\n"
          "as a function of *f_grid*, *p_grid*, *lat_grid*, and *lon_grid*. \n"
          "\n"
          "This is mainly for testing and plotting gas absorption. For RT\n"
          "calculations, gas absorption is calculated or extracted locally,\n"
          "therefore there is no need to store a global field. But this variable\n"
          "is handy for easy plotting of absorption vs. pressure, for example.\n"
          "\n"
          "Unit:       1/m\n"
          "\n"
          "Dimensions: [species, f_grid, *stokes_dim*, stokes_dim, p_grid, lat_grid, lon_grid]\n"),
      GROUP("Tensor7")));

  wsv_data.push_back(WsvRecord(
      NAME("psd_data"),
      DESCRIPTION(
          "Particle size distribution values for a set of points.\n"
          "\n"
          "The variable contains particle size distribution data for one scattering\n"
          "species. The row dimension corresponds to different positions, in the\n"
          "same way as *pnd_agenda_input* is defined.\n"
          "\n"
          "Dimensions: [ n_points, n_scattering_elements ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("psd_size_grid"),
      DESCRIPTION(
          "The particle sizes associated with *psd_data*.\n"
          "\n"
          "This variable holds the size of each scattering element considered.\n"
          "Size can be defined differently, depending on particle size distribution\n"
          "used. Most common choices should by equivalent diameter, maximum diameter\n"
          "and mass.\n"
          "\n"
          "Dimension: [ n_sizes ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("p_grid"),
      DESCRIPTION(
          "The pressure grid.\n"
          "\n"
          "The pressure levels on which the atmospheric fields are defined.\n"
          "This variable must always be defined. The grid must be sorted in\n"
          "decreasing order, with no repetitions.\n"
          "\n"
          "No gap between the lowermost pressure level and the surface is \n"
          "allowed. The uppermost pressure level defines the practical upper\n"
          "limit of the atmosphere as vacuum is assumed above.\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  Pa\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("p_grid_orig"),
      DESCRIPTION(
          "The original pressure grid before optimization.\n"
          "\n"
          "This variable is used to interpolate *cloudbox_field* back to its original\n"
          "size after the calculation with *OptimizeDoitPressureGrid*.\n"
          "The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  Pa\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("p_hse"),
      DESCRIPTION(
          "Reference pressure calculation of hydrostatic equilibrium.\n"
          "\n"
          "The altitude specified by this pressure is used as the reference\n"
          "when calculating hydrostatic equilibrium. That is, the geometrical\n"
          "altitude at this pressure is not changed.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  Pa\n"),
      GROUP("Numeric")));
  wsv_data.push_back(WsvRecord(
      NAME("radiance_field"),
      DESCRIPTION(
          "Radiance field.\n"
          "\n"
          "Radiant flux received by a surface per unit solid angle and per unit\n"
          "area for each hemisphere. The last dimension denotes the hemispheres.\n"
          "The first component is the downward radiance and the second component\n"
          "is the upward radiance.\n"
          "\n"
          "Units: W / (m^2 sr)\n"
          "\n"
          "Size: [p_grid, \n"
          "       lat_grid, \n"
          "       lon_grid, \n"
          "       N_za, N_aa\n"),
      GROUP("Tensor5")));

  wsv_data.push_back(WsvRecord(
      NAME("range_bins"),
      DESCRIPTION(
          "The range bins of an active instrument.\n"
          "\n"
          "The bins are assumed to cover a range without gaps, and the bins are\n"
          "defined by their edges. That is, the length of this vector is the\n"
          "number of bins + 1.\n"
          "\n"
          "The bins can potentially be defined in two ways, by altitude or time.\n"
          "See the method you are using, if this variable shall hold time or\n"
          "altitude (or maybe both options are treated).\n"
          "\n"
          "Unit: m or s\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("refr_index_air"),
      DESCRIPTION(
          "Real part of the refractive index of air.\n"
          "\n"
          "The variable contains the refractive index summed over all relevant\n"
          "constituents, at one position in the atmosphere. This refractive\n"
          "is related to the phase velocity. See also *refr_index_air_group*.\n"
          "\n"
          "Unit: 1\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("refr_index_air_agenda"),
      DESCRIPTION("Agenda calculating the refractive index of air.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("refr_index_air_group"),
      DESCRIPTION(
          "Group index of refractivity.\n"
          "\n"
          "This variable is defined as the ratio between group velocity and the\n"
          "speed of ligh in vacuum. That is, it is defined as the \"standard\"\n"
          "refractive index, but refers to the group velocity instead of the\n"
          "phase velocity. See also *refr_index_air*.\n"
          "\n"
          "Unit: 1\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("refellipsoid"),
      DESCRIPTION(
          "Reference ellipsoid.\n"
          "\n"
          "This vector specifies the shape of the reference ellipsoid. The\n"
          "vector must have length 2, where the two elements are:\n"
          "\n"
          "1. Equatorial radius.\n"
          "2. The eccentricity.\n"
          "\n"
          "The eccentricity is sqrt(1-b * b/a * a) where a and b are equatorial and\n"
          "polar radius, respectively. If the eccentricity is set to 0, an\n"
          "average radius should be used instead of the equatorial one.\n"
          "\n"
          "The eccentricity must be 0 for 1D calculations, as a spherical Earth\n"
          "is implied by setting *atmosphere_dim* to 1. For 2D, the selected\n"
          "ellipsoid parameters should be selected according to cross-section\n"
          "between the real ellipsoid and the 2D plane considered. That is\n"
          "the applied ellipsoid shall have een converted to match the internal\n"
          "treatment of 2D cases. For 3D, models can be used, such as WGS84.\n"
          "\n"
          "Usage:  Set by the user.\n"
          "\n"
          "Size:   [ 2 ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("retrieval_checked"),
      DESCRIPTION(
          "Flag indicating completeness and consistency of retrieval setup.\n"
          "\n"
          "Unit: Boolean\n"),
      GROUP("Index"), Index{0}));

  wsv_data.push_back(WsvRecord(
      NAME("retrieval_eo"),
      DESCRIPTION(
          "The estimated error in the retrieval due to uncertainty in the observations.\n"
          "\n"
          "The vector contains the square roots  of the diagonal elements of  the\n"
          "covariance matrix of the error due to measurement noise, S_m in Rodgers'\n"
          "book.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("retrieval_ss"),
      DESCRIPTION(
          "The estimated error in the retrieval due to limited resolution of the\n"
          "observation system.\n"
          "\n"
          "The vector contains the square roots of the diagonal\n"
          "elements of the covariance matrix of the smoothing error, S_s in Rodgers'\n"
          "book.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("rte_alonglos_v"),
      DESCRIPTION(
          "Velocity along the line-of-sight to consider for a RT calculation.\n"
          "\n"
          "This variable gives the velocity of the imaginary detector in\n"
          "monochromatic pencil beam calculations. The relevant velocity is\n"
          "the projection along the line-of-sight (ie. total velocity shall not\n"
          "be given). A positive value means a movement of the detector in the\n"
          "same direction as the line-of-sight.\n"
          "\n"
          "This variable is required to include Doppler effects due to\n"
          "velocities of the observer, relative the centre of the coordinate\n"
          "system used that is fixed to the planets centre point.\n"
          "\n"
          "Unit: [ m/s ]\n"),
      GROUP("Numeric"), Numeric{0.0}));

  wsv_data.push_back(WsvRecord(
      NAME("rte_los"),
      DESCRIPTION(
          "A line-of-sight for (complete) radiative transfer calculations.\n"
          "\n"
          "This variable gives the observation direction for monochromatic\n"
          "pencil beam calculations. Hence, it is the line-of-sight at the end\n"
          "point of the propagation path.\n"
          "\n"
          "For 1D and 2D cases, *rte_los* is a vector of length 1 holding the \n"
          "zenith angle. For 3D, the length of the vector is 2, where the\n"
          "additional element is the azimuthal angle. These angles are defined\n"
          "in the ARTS user guide (AUG). Look in the index for \"zenith angle\"\n"
          "and \"azimuthal angle\".\n"
          "\n"
          "Usage: See above.\n"
          "\n"
          "Units: [ degree, degree ]\n"
          "\n"
          "Size:  [ 1 or 2 ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("rte_pos"),
      DESCRIPTION(
          "A geographical position for starting radiative transfer calculations.\n"
          "\n"
          "This variable gives the observation position for monochromatic\n"
          "pencil beam calculations. Hence, it is the end point of the\n"
          "propagation path.\n"
          "\n"
          "This variable is a vector with a length equalling the atmospheric\n"
          "dimensionality. The first element is the geometrical altitude.\n"
          "Element 2 is the latitude and element 3 is the longitude.\n"
          "\n"
          "Usage: See above. \n"
          "\n"
          "Units: [ m, degree, degree ]\n"
          "\n"
          "Size:  [ atmosphere_dim ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("rte_pos2"),
      DESCRIPTION(
          "A second geographical position to define the geometry for\n"
          "radiative transfer calculations.\n"
          "\n"
          "This variable is used when the propagation path is defined by two\n"
          "positions, instead of a position (*rte_pos*) and a line-of-sight\n"
          "(*rte_los*). That is, this variable basically replaces *rte_los*\n"
          "for the cases of consideration. In practice, *rte_los* is determined\n"
          "by finding the propagation path between *rte_pos* and *rte_pos2*.\n"
          "\n"
          "As *rte_pos* with the exception that a \"latitude\" must also be\n"
          "specified for 1D. This is the angular distance to *rte_pos*, where\n"
          "this distance is defined as the 2D-\"latitude\".\n"
          "\n"
          "Usage: See above. \n"
          "\n"
          "Units: [ m, degree, degree ]\n"
          "\n"
          "Size:  [ atmosphere_dim ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("rtp_mag"),
      DESCRIPTION(
          "Magnetic field at a radiative transfer point.\n"
          "\n"
          "See *mag_u_field* etc. for a definition of the different components.\n"
          "For this variable the components are put together and thus defines\n"
          "magnetic field vector. Hence, this is a vector of length three, even\n"
          "if any of the input fields is set to be empty.\n"
          "\n"
          "The WSV is used as input to methods and agendas calculating radiative\n"
          "properties for a given conditions.\n"
          "\n"
          "Usage: Communication variable.\n"
          "\n"
          "Units: T\n"
          "\n"
          "Size:  [ u-component, v-component, w-component ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("rtp_los"),
      DESCRIPTION(
          "Line-of-sight at a radiative transfer point.\n"
          "\n"
          "This variable holds a local line-of-sight. The angles of this\n"
          "vector are defined as for *rte_los*.\n"
          "\n"
          "The WSV is used as input to methods and agendas calculating radiative\n"
          "properties for a given conditions.\n"
          "\n"
          "Usage: Communication variable.\n"
          "\n"
          "Units: [ degree, degree ]\n"
          "\n"
          "Size:  [ 1 or 2 ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("rtp_pos"),
      DESCRIPTION(
          "Position of a radiative transfer point.\n"
          "\n"
          "This vector is defined as *rte_pos*, but holds a position along\n"
          "the propgation path, or the start point for new paths, in contrast\n"
          "to *rte_pos* that is position of the (imaginary) detector.\n"
          "\n"
          "The WSV is used as input to methods and agendas calculating radiative\n"
          "properties for a given conditions.\n"
          "\n"
          "Usage: Communication variable.\n"
          "\n"
          "Units: [ m, degree, degree ]\n"
          "\n"
          "Size:  [ atmosphere_dim ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("rtp_pressure"),
      DESCRIPTION(
          "Pressure at a radiative transfer point.\n"
          "\n"
          "This scalar variable holds the local pressure.\n"
          "\n"
          "The WSV is used as input to methods and agendas calculating radiative\n"
          "properties for a given conditions.\n"
          "\n"
          "Usage: Communication variable.\n"
          "\n"
          "Units: [ Pa ]\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("rtp_temperature"),
      DESCRIPTION(
          "Temperature at a radiative transfer point.\n"
          "\n"
          "This scalar variable can hold the local temperature. It is intended\n"
          "mainly for communication with various methods and agendas, such as\n"
          "methods and agendas calculating absorption coefficients.\n"
          "The WSV is used as input to methods and agendas calculating radiative\n"
          "properties for a given conditions.\n"
          "\n"
          "Usage: Communication variable.\n"
          "\n"
          "Units: [ K ]\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("rt_integration_option"),
      DESCRIPTION(
          "Switch between integration approaches for radiative transfer steps.\n"
          "\n"
          "See each WSM using this varaible as input for available options.\n"),
      GROUP("String"), String{"default"}));

  wsv_data.push_back(WsvRecord(
      NAME("rtp_nlte"),
      DESCRIPTION(
          "NLTE temperature/ratio at a radiative transfer point.\n"
          "\n"
          "This vector variable can hold the NLTE temperature/ratio. It is intended\n"
          "mainly for communication with various methods and agendas, such as\n"
          "methods and agendas calculating absorption coefficients.\n"
          "The WSV is used as input to methods and agendas calculating radiative\n"
          "properties for a given conditions.\n"
          "\n"
          "Usage: Communication variable.\n"
          "\n"
          "Units: [ K/# ]\n\n"
          "Size:  [ NLTE levels, 1, 1, 1 ] or [ 0, 0, 0, 0 ]\n"),
      GROUP("EnergyLevelMap")));

  wsv_data.push_back(WsvRecord(
      NAME("rtp_vmr"),
      DESCRIPTION(
          "Absorption species abundances for radiative transfer calculations.\n"
          "\n"
          "This vector variable holds the local abundance of the constituents\n"
          "included in *abs_species*.\n"
          "\n"
          "The WSV is used as input to methods and agendas calculating radiative\n"
          "properties for a given conditions.\n"
          "\n"
          "Usage: Communication variable.\n"
          "\n"
          "Units: [ Differ between the elements, can be VMR, kg/m3 or #/m3. ]\n"
          "\n"
          "Size:  Should match abs_species.nelem()\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("scat_data"),
      DESCRIPTION(
          "Array of single scattering data.\n"
          "\n"
          "As *scat_data_raw*, but with frequency grids and dimensions reduced\n"
          "to the RT's *f_grid* or a single frequency entry. Also, temperature\n"
          "grid or dimensions can be reduced to a single entry, meaning no\n"
          "temperature interpolation is done for the respective data.\n"
          "\n"
          "Standard approach to derive scat_data is to use *scat_dataCalc* to\n"
          "derive it from *scat_data_raw*."),
      GROUP("ArrayOfArrayOfSingleScatteringData")));

  wsv_data.push_back(WsvRecord(
      NAME("scat_data_checked"),
      DESCRIPTION(
          "OK-flag for *scat_data*.\n"
          "\n"
          "Relevant checks are performed by *scat_data_checkedCalc*. Only the\n"
          "value 1 is taken as OK.\n"),
      GROUP("Index"), Index{0}));

  wsv_data.push_back(WsvRecord(
      NAME("scat_data_raw"),
      DESCRIPTION(
          "Array of raw single scattering data.\n"
          "\n"
          "This variable holds the single scattering properties for all \n"
          "scattering elements, organized according to their assignment to a\n"
          "scattering species. *scat_data_raw* entries can be derived from\n"
          "precalculated data files using the methods *ScatElementsPndAndScatAdd*,\n"
          "*ScatSpeciesPndAndScatAdd*, or *ScatSpeciesScatAndMetaRead* or\n"
          "can be calculated using *scat_data_singleTmatrix*.\n"
          "\n"
          "This may be used in combination with *scat_meta*\n"
          "\n"
          "Usage: Method ouput.\n"
          "\n"
          "Members:\n"
          "\n"
          "- SingleScatteringData:\n"
          "\n"
          "  - Enum[ptype attribute]\n"
          "  - String[description] \n"
          "  - Vector[f_grid]\n"
          "  - Vector[T_grid]\n"
          "  - Vector[za_grid]\n"
          "  - Vector[aa_grid]\n"
          "  - Tensor7[pha_mat_data]\n"
          "\n"
          "    - [f_grid, T_grid, za_grid(scattered), aa_grid(scattered), za_grid(incoming), aa_grid(incoming), matrix_element]\n"
          "\n"
          "  - Tensor5[ext_mat_data]\n"
          "\n"
          "    - [f_grid, T_grid, za_grid, aa_grid, matrix_element]\n"
          "\n"
          "  - Tensor5[abs_vec_data]\n"
          "\n"
          "    - [f_grid, T_grid, za_grid, aa_grid, matrix_element]\n"
          "\n"
          "Dimensions: [number of scattering species][number of scattering elements] \n"),
      GROUP("ArrayOfArrayOfSingleScatteringData")));

  wsv_data.push_back(WsvRecord(
      NAME("scat_data_mono"),
      DESCRIPTION(
          "Monochromatic single scattering data.\n"
          "\n"
          "This variable holds the single scattering properties for all\n"
          "scattering species and scattering elements for a specified frequency.\n"
          "It can be calculated from *scat_data* using *scat_data_monoCalc*,\n"
          "which interpolates *scat_data* to the required frequency.\n"),
      GROUP("ArrayOfArrayOfSingleScatteringData")));

  wsv_data.push_back(WsvRecord(
      NAME("scat_data_single"),
      DESCRIPTION(
          "Structure for the single scattering data.\n"
          "\n"
          "Comprises the single scattering data of a single scattering element.\n"
          "See ARTS user guide for further information.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Dimensions:\n"
          "\n"
          "- SingleScatteringData \n"
          "\n"
          "  - Enum[ptype attribute]\n"
          "  - String[description] \n"
          "  - Vector[f_grid]\n"
          "  - Vector[T_grid]\n"
          "  - Vector[za_grid]\n"
          "  - Vector[aa_grid]\n"
          "  - Tensor7[pha_mat_data]\n"
          "\n"
          "    - [f_grid, T_grid, za_grid(scattered), aa_grid(scattered), za_grid(incoming), aa_grid(incoming), matrix_element]\n"
          "\n"
          "  - Tensor5[ext_mat_data]\n"
          "\n"
          "    - [f_grid, T_grid, za_grid, aa_grid, matrix_element]\n"
          "\n"
          "  - Tensor5[abs_vec_data]\n"
          "\n"
          "    - [f_grid, T_grid, za_grid, aa_grid, matrix_element]\n"),
      GROUP("SingleScatteringData")));

  wsv_data.push_back(WsvRecord(
      NAME("scat_lat_index"),
      DESCRIPTION(
          "Latitude index for scattering calculations.\n"
          "\n"
          "This variable is used in methods used for computing scattering\n"
          "properties of scattering elements like *opt_prop_sptFromData* and\n"
          "*pha_matCalc*. It holds the information about the position for which the\n"
          "scattering calculations are done.\n"
          "\n"
          "Usage:    Input to the methods *spt_calc_agenda*, *pha_mat_spt_agenda*\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("scat_lon_index"),
      DESCRIPTION(
          "Longitude index for scattering calculations.\n"
          "\n"
          "This variable is used in methods used for computing scattering\n"
          "properties of scattering elements like *opt_prop_sptFromData* and\n"
          "*pha_matCalc*. It holds the information about the position for which the\n"
          "scattering calculations are done.\n"
          "\n"
          "Usage:    Input to the methods *spt_calc_agenda*, *pha_mat_spt_agenda*\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("scat_meta_single"),
      DESCRIPTION(
          "Structure for the scattering meta data.\n"
          "\n"
          "This variable holds the scattering meta data for a single scattering\n"
          "element (see AUG for definition). Scattering meta data comprises\n"
          "the microphysical description of the scattering element as necessary\n"
          "to relate single scattering properties with mass density or flux\n"
          "fields. That is, e.g., in order to handle the scattering element in\n"
          "particle size (and shape) distribution calculations.\n"
          "\n"
          "For a definition of the structure members see below.\n"
          "\n"
          "Members of Numeric type can be flagged as unknown by setting them to\n"
          "NAN. This will cause a runtime error in case the parameter is needed in\n"
          "the calculation, but will be ignored otherwise.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Members:\n"
          "  description [*String*]\n"
          "    Description: Free-form description of the scattering element,\n"
          "    holding information deemed of interest by the user but not covered\n"
          "    by other structure members (and not used within ARTS).\n"
          "  source [*String*]\n"
          "    Description: Free-form description of the source of the data,\n"
          "    e.g., Mie, T-Matrix, or DDA calculation or a database or a\n"
          "    literature source.\n"
          "  refr_index [*String*]\n"
          "    Description: Free-form description of the underlying complex\n"
          "    refractive index data, e.g., a literature source.\n"
          "  mass [*Numeric*]\n"
          "    Unit: [kg]\n"
          "    Description: The mass of the scattering element.\n"
          "  diameter_max [*Numeric*]\n"
          "    Unit: [m]\n"
          "    Description: The maximum diameter (or dimension) of the scattering\n"
          "    element, defined by the circumferential sphere diameter of the\n"
          "    element. Note that this parameter is only used by some size\n"
          "    distributions; it does not have a proper meaning if the scattering\n"
          "    element represents an ensemble of differently sized particles.\n"
          "  diameter_volume_equ [*Numeric*]\n"
          "    Unit: [m]\n"
          "    Description: The volume equivalent sphere diameter of the\n"
          "    scattering element, i.e., the diameter of a sphere with the same\n"
          "    volume. For nonspherical particles, volume refers to the volume\n"
          "    of the particle-forming substance, not that of the circumferential\n"
          "    sphere (which can be derived from diameter_max). If the particle\n"
          "    consists of a mixture of materials, the substance\n"
          "    encompasses the complete mixture. E.g., the substance of 'soft'\n"
          "    ice particles includes both the ice and the air.\n"
          "  diameter_area_equ_aerodynamical [*Numeric*]\n"
          "    Unit: [m]\n"
          "    Description: The area equivalent sphere diameter of the\n"
          "    scattering element, i.e., the diameter of a sphere with the same\n"
          "    cross-sectional area. Here, area refers to the aerodynamically\n"
          "    relevant area, i.e., the cross-sectional area perpendicular to the\n"
          "    direction of fall. Similarly to volume in the definition of\n"
          "    diameter_volume_equ, for non-spherical and mixed-material\n"
          "    particles, area refers to the area covered by the substance\n"
          "    mixture of the particle.\n"),
      GROUP("ScatteringMetaData")));

  wsv_data.push_back(WsvRecord(
      NAME("scat_meta"),
      DESCRIPTION(
          "An Array of scattering meta data (*scat_meta_single*).\n"
          "\n"
          "The array holds the meta data for all scattering elements. For a\n"
          "description of the meta data contents refer to the documentation\n"
          "of *scat_data_single*.\n"
          "\n"
          "Corresponding to *scat_data*, it is organized in terms of scattering\n"
          "species (i.e., one sub-array per scattering species holding one\n"
          "*scat_meta_single* instance per scattering element assigned to this\n"
          "scattering species). It is primarily used for particle size and shape\n"
          "distribution calculations using *pnd_fieldCalcFromParticleBulkProps*.\n"
          "It is also applied for deducing microphysical characterizations of\n"
          "scattering species, e.g., by *particle_massesFromMetaData*.\n"
          "\n"
          "Note: This array must contain as many elements as *scat_data* (on\n"
          "both array levels).\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Dimensions: [scattering species][scattering elements]\n"
          "\n"
          "See also: *scat_meta_single*.\n"),
      GROUP("ArrayOfArrayOfScatteringMetaData")));

  wsv_data.push_back(WsvRecord(
      NAME("scat_p_index"),
      DESCRIPTION(
          "Pressure index for scattering calculations.\n"
          "\n"
          "This variable is used in methods used for computing scattering\n"
          "properties of scattering elements like *opt_prop_sptFromData* and\n"
          "*pha_matCalc*. It holds the information about the location for which the\n"
          "scattering calculations are done.\n"
          "\n"
          "Usage:    Input to the methods *spt_calc_agenda*, *pha_mat_spt_agenda*\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("scat_species"),
      DESCRIPTION(
          "Array of Strings defining the scattering species to consider.\n"
          "\n"
          "Each String contains the information to connect scattering species\n"
          "(e.g., hydrometeor) atmospheric fields with the microphysical\n"
          "information like size and shape distributions. The strings follow\n"
          "the following structure with individual elements separated by dashes:\n"
          "\n"
          "- scattering species name [*String*]\n"
          "  the name of the scattering species' atmospheric field. Free form,\n"
          "  but is matched to *atm_fields_compact* fields by their names.\n"
          "  Common are, e.g., IWC (ice water content), LWC (liquid water\n"
          "  content), RR (rain rate), and SR (snow rate).\n"
          "- particle size distribution [*String*]:\n"
          "  the size distribution function/parametrization to apply. For\n"
          "  currently possible PSDs see *pnd_fieldCalcFromParticleBulkProps*.\n"
          "\n"
          "Example: [''IWC-MH97'', ''LWC-H98_STCO'', ...]\n"),
      GROUP("ArrayOfString"),
      ArrayOfString{}));

  wsv_data.push_back(WsvRecord(
      NAME("scat_species_a"),
      DESCRIPTION(
          "Mass-size relationship parameter, for one scattering species.\n"
          "\n"
          "Some methods require a relationship between mass and particle size,\n"
          "valid for the complete scattering species. A common model for this\n"
          "relationship is::\n"
          "\n"
          "  mass(x) = a * x^b,\n"
          "\n"
          "where x is size (that could be Dveq, Dmax or mass) and a/b are parameters.\n"
          "\n"
          "This WSV is a in the expression above.\n"
          "The WSV matching b is *scat_species_b*.\n"
          "The WSV matching x is *scat_species_x*.\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("scat_species_b"),
      DESCRIPTION(
          "Mass-size relationship parameter, for one scattering species.\n"
          "\n"
          "See *scat_species_a* for details.\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("scat_species_x"),
      DESCRIPTION(
          "The size grid of one scattering species.\n"
          "\n"
          "The variable holds the sizes associated with one scattering species.\n"
          "The typical application of these data are as the size grid when\n"
          "calculating particle size distributions.\n"
          "\n"
          "The user must set this WSV as several quantities can be used as size,\n"
          "such as mass and maximum diamater.\n"
          "\n"
          "See also *scat_species_a*, for example usage of this WSV.\n"
          "\n"
          "Dimension:  [number of scattering elements]\n"),
      GROUP("Vector")));

  wsv_data.push_back(
      WsvRecord(NAME("select_abs_species"),
                DESCRIPTION(R"--(A select species tag group from *abs_species*

If set to empty, this selection is void.  It must otherwise match perfectly a tag inside
*abs_species* for that to be the selection.
)--"),
                GROUP("ArrayOfSpeciesTag"), ArrayOfSpeciesTag{}));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_checked"),
      DESCRIPTION(
          "OK-flag for sensor related variables.\n"
          "\n"
          "This variable flags that sensor variables are defined in a formally\n"
          "and practically correct way. For example, it checks for correct\n"
          "dimensions of *sensor_pos* and *sensor_los*.\n"
          "\n"
          "Shall be set by *sensor_checkedCalc*. See that WSM for treated WSVs.\n"
          "Only the value 1 is taken as OK.\n"),
      GROUP("Index"), Index{0}));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_description_amsu"),
      DESCRIPTION(
          "Sensor description for simple AMSU setup.\n"
          "\n"
          "This is a compact description of an AMSU-type sensor. The matrix\n"
          "contains one row for each instrument channel. Each row contains three\n"
          "elements: LO position [Hz], offset of the channel center from the LO\n"
          "[Hz], and channel width [Hz].\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit: All entries in Hz.\n"
          "\n"
          "Size: [number of channels, 3]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_los"),
      DESCRIPTION(
          "The sensor line-of-sight (LOS) for each measurement block.\n"
          "\n"
          "Line-of-sights are specified by giving the zenith and azimuth angles.\n"
          "Column 1 holds the zenith angle. This angle is simply the angle \n"
          "between the zenith and LOS directions. For 1D and 3D the valid\n"
          "range is [0 180], while for 2D angles down to -180 degrees are\n"
          "allowed. Negative angles signifies for 2D observations towards\n"
          "lower latitudes, while positive angles means observations towards\n"
          "higher latitudes. Nadir corresponds throughout to 180 degrees.\n"
          "\n"
          "The azimuth angle is given with respect to the meridian plane. That\n"
          "is, the plane going through the north and south poles. The valid \n"
          "range is [-180,180] where angles are counted clockwise; 0 means\n"
          "that the viewing or propagation direction is north-wise and +90 means\n"
          "that the direction of concern goes eastward.\n"
          "\n"
          "No azimuth angle shall be specified for 1D and 2D. This angle is in\n"
          "general of no concern for these atmospheric dimensionalities, but\n"
          "matter in some cases, such as with respect to the Doppler shift due\n"
          "to winds. For 1D the azimuth angle is then assumed to be 0 deg, i.e.\n"
          "the sensor is treated to be directed towards North. For 2D, the \n"
          "implied azimuth is 0 or 180, depending of the zenith angle is positive\n"
          "or negative.\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  [ degrees, degrees ]\n"
          "\n"
          "Size:  [ number of measurement blocks, 1 or 2 ]\n"),
      GROUP("Matrix")));

    wsv_data.push_back(WsvRecord(
      NAME("sensor_los_geodetic"),
      DESCRIPTION(
          "As *sensor_los* but matching geodetic coordinates.\n"
          "\n"
          "For this version zenith is defined as the normal of the reference\n"
          "ellipsoid, in contrast to *sensor_los* zenith is along the direction\n"
          "towards the planets centre.\n" 
          "\n"
          "Probably only useful for 3D.\n"
          "\n"          
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  [ degrees, degrees ]\n"
          "\n"
          "Size:  [ number of measurement blocks, 2 ]\n"),
      GROUP("Matrix")));

    wsv_data.push_back(WsvRecord(
      NAME("sensor_los_ecef"),
      DESCRIPTION(
          "As *sensor_los* but matching ECEF coordinates.\n"
          "\n"
          "For this version of sensor_los, each row shall hold [dx,dy,dz],\n"
          "where dx, dy and dz are the x, y and z components if the line-of-sight\n"
          "directions in ECEF coordinates. [dx,dy,dz] must form a unit vector (i.e.\n" 
          "its 2-norm shall be 1).\n"
          "\n"
          "Probably only useful for 3D.\n"
          "\n"          
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  [ m, m, m ]\n"
          "\n"
          "Size:  [ number of measurement blocks, 3 ]\n"),
      GROUP("Matrix")));
    
  wsv_data.push_back(WsvRecord(
      NAME("sensor_norm"),
      DESCRIPTION(
          "Flag if sensor response should be normalised or not (0 or 1).\n"
          "\n"
          "If the flag is set to 1 each sensor response is normalised (where\n"
          "applicable). If set to 0 the sensor responses are left as provided.\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a sub-entry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage: Set by the user.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_pol"),
      DESCRIPTION(
          "A set of polarisation response angles.\n"
          "\n"
          "The standard choice to consider the polarisation response of the\n"
          "reciever is by *instrument_pol*, and this response becomes then part\n"
          "of *sensor_response*. However, that choice is not possible when the\n"
          "polartisation response changes between measurement blocks, and this\n"
          "variable combined with the *yApplySensorPol* offers an alternative for\n"
          "such situations. This WSV also allows defintion of an arbitrary\n"
          "polarisation angle.\n"
          "\n"
          "When applying the polarisation response by *yApplySensorPol*, this\n"
          "variable complements *sensor_pos* and *sensor_los*. This WSV matrix\n"
          "is also a matrix, that shall have the same number of rows as the other\n"
          "two matrices. \n"
          "\n"
          "The columns of *sensor_pol* corresponds to the channels/frequencies\n"
          "of the receiver. Each element gives the polarisation angle. A pure\n"
          "vertical response has the angle 0 deg, and pure horisontal 90 deg.\n"
          "If all U values (Stokes element 3) are zero, the sign of the angle does,\n"
          "not matter, and 0 and 180 degrees give the same result. With non-zero\n"
          "U, the result of e.g. -45 and +45 degrees differ.\n"
          "\n"
          "Note that a receiver with a linear response is assumed. Circular\n"
          "polarisation is not affected by any rotation.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  [ degrees ]\n"
          "\n"
          "Size:  [ number of measurement blocks, number of channels/frequencies ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_pos"),
      DESCRIPTION(
          "The sensor position for each measurement block.\n"
          "\n"
          "The sensor positions are specified as a matrix, where the number of\n"
          "columns shall be equal to *atmosphere_dim*. Column 1 shall contain\n"
          "the altitude of the sensor platform, column 2 the latitude and the \n"
          "last column the longitude. The number of rows corresponds to the\n"
          "number of measurement blocks.\n"
          "\n"
          "Valid range for latitudes in 3D is [-90,90], while for 2D any value\n"
          "is accepted. Accepted range for longitudes are [-360,360].\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  [ m, degrees, degrees ]\n"
          "\n"
          "Size:  [ number of measurement blocks, atmosphere_dim ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_pos_geodetic"),
      DESCRIPTION(
          "As *sensor_pos* but using geodetic coordinates.\n"
          "\n"
          "For this version the second column shall hold geodetic latitudes,\n"
          "in contrast to *sensor_pos* where the geocentric system us used.\n"
          "Please note that also the altitude (column 1) differs between\n"
          "the two versions of the variables. Here the altitude is with\n"
          "taken along the local nadir, while for *sensor_pos* it is taken\n"  
          "along the direction towards the planets centre.\n" 
          "\n"
          "Probably only useful for 3D.\n"
          "\n"          
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  [ m, degrees, degrees ]\n"
          "\n"
          "Size:  [ number of measurement blocks, atmosphere_dim ]\n"),
      GROUP("Matrix")));
  
  wsv_data.push_back(WsvRecord(
      NAME("sensor_pos_ecef"),
      DESCRIPTION(
          "As *sensor_pos* but using ECEF coordinates.\n"
          "\n"
          "The sensor position is here specified as earth-centered, earth-fixed\n"
          "(ECEF) coordinates (using standard definition of ECEF).\n"
          "\n"
          "Probably only useful for 3D.\n"
          "\n"
          "Column 1, 2 and 3 shall hold x, y and z coordinate, respectively.\n"
          "\n"          
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  [ m, m, m ]\n"
          "\n"
          "Size:  [ number of measurement blocks, atmosphere_dim ]\n"),
      GROUP("Matrix")));
  
  wsv_data.push_back(WsvRecord(
      NAME("sensor_response"),
      DESCRIPTION(
          "The matrix modelling the total sensor response.\n"
          "\n"
          "This matrix describes the sensor respons for one measurement block\n"
          "The response is assumed to be identical for each such block.\n"
          "\n"
          "The matrix is the product of all the individual sensor response\n"
          "matrices. Therefore its dimensions are depending on the total sensor\n"
          "configuration. The *sensor_response* has to initialised by the \n"
          "*sensor_responseInit* method.\n"
          "\n"
          "Usage: Output/input to the ``sensor_response...`` methods.\n"
          "\n"
          "Units: -\n"
          "\n"
          "Dimension: See the individual ``sensor_response...`` method. \n"),
      GROUP("Sparse")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_response_aa"),
      DESCRIPTION(
          "The relative azimuth angles associated with the output of\n"
          "*sensor_response*.\n"
          "\n"
          "The variable shall not be set manually, it will be set together with\n"
          "*sensor_response* by sensor response WSMs.\n"
          "\n"
          "Usage: Set by sensor response methods.\n"
          "\n"
          "Unit:  [ degrees ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_response_agenda"),
      DESCRIPTION(
          "Agenda providing the sensor response data for a measurement block.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_response_dlos"),
      DESCRIPTION(
          "The relative zenith and azimuth angles associated with the output of\n"
          "*sensor_response*.\n"
          "\n"
          "Definition of angles match *mblock_dlos*. Works otherwise as\n"
          "*sensor_response_f*.\n"
          "\n"
          "The variable shall not be set manually, it will be set together with\n"
          "*sensor_response* by sensor response WSMs.\n"
          "\n"
          "Usage: Set by sensor response methods.\n"
          "\n"
          "Unit:  [ degrees ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_response_dlos_grid"),
      DESCRIPTION(
          "The zenith and azimuth angles associated with *sensor_response*.\n"
          "\n"
          "A variable for communication between sensor response WSMs. Matches\n"
          "initially *mblock_dlos*, but is later adjusted according to the\n"
          "sensor specifications. Only defined when a common grid exists. Values\n"
          "are here not repeated as in *sensor_response_dlos*.\n"
          "\n"
          "Usage: Set by sensor response methods.\n"
          "\n"
          "Unit:  [ degrees ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_response_f"),
      DESCRIPTION(
          "The frequencies associated with the output of *sensor_response*.\n"
          "\n"
          "This vector gives the frequency for each element of the measurement\n"
          "vector produced inside one measurement block. The frequencies of\n"
          "the total measurement vector, *y*, are obtained by repeating these\n"
          "frequencies n times, where n is the number of measurement blocks\n"
          "(e.g. the number of rows in *sensor_pos*).\n"
          "\n"
          "The variable shall not be set manually, it will be set together with\n"
          "*sensor_response* by sensor response WSMs.\n"
          "\n"
          "Usage: Set by sensor response methods.\n"
          "\n"
          "Unit:  [ Hz ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_response_f_grid"),
      DESCRIPTION(
          "The frequency grid associated with *sensor_response*.\n"
          "\n"
          "A variable for communication between sensor response WSMs. Matches\n"
          "initially *f_grid*, but is later adjusted according to the sensor\n"
          "specifications. Only defined when a common grid exists. Values are\n"
          "here not repeated as in *sensor_response_f*\n"
          "\n"
          "Usage: Set by sensor response methods.\n"
          "\n"
          "Unit:  [ Hz ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_response_pol"),
      DESCRIPTION(
          "The polarisation states associated with the output of\n"
          "*sensor_response*.\n"
          "\n"
          "Works basically as *sensor_response_f*.\n"
          "\n"
          "See *instrument_pol* for coding of polarisation states.\n"
          "\n"
          "The variable shall not be set manually, it will be set together with\n"
          "*sensor_response* by sensor response WSMs.\n"
          "\n"
          "Usage: Set by sensor response methods.\n"
          "\n"
          "Unit:  [ - ]\n"),
      GROUP("ArrayOfIndex")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_response_pol_grid"),
      DESCRIPTION(
          "The \"polarisation grid\" associated with *sensor_response*.\n"
          "\n"
          "A variable for communication between sensor response WSMs. It is\n"
          "initially 1:stokes_dim, but can later adjusted according to the \n"
          "sensor specifications. Only defined when a common grid exists. \n"
          "\n"
          "See *instrument_pol* for coding of polarisation states.\n"
          "\n"
          "Usage: Set by sensor response methods.\n"
          "\n"
          "Unit:  [ - ]\n"),
      GROUP("ArrayOfIndex")));

  wsv_data.push_back(WsvRecord(
      NAME("sensor_time"),
      DESCRIPTION(
          "The time for each measurement block.\n"
          "\n"
          "This WSV is used when a time must be assigned to the measurements.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  [ UTC date and time ]\n"
          "\n"
          "Size:  [ number of measurement blocks ]\n"),
      GROUP("ArrayOfTime")));

  wsv_data.push_back(WsvRecord(
      NAME("sideband_mode"),
      DESCRIPTION(
          "Description of target sideband.\n"
          "\n"
          "A text string describing which of the two sidebands (of a heterodyne\n"
          "instrument) that can be seen as \"main\" band. Possible choices are:\n"
          "\n"
          "- \"lower\" : Low frequency sideband shall be considered as target.\n"
          "- \"upper\" : High frequency sideband shall be considered as target.\n"
          "\n"
          "Usage: Set by the user.\n"),
      GROUP("String")));

  wsv_data.push_back(WsvRecord(
      NAME("sideband_mode_multi"),
      DESCRIPTION(
          "Description of target sideband for a multiple LO receiver.\n"
          "\n"
          "As *sideband_mode* but handles an instrument with several LO chains.\n"
          "See further *lo_multi* and *sideband_response_multi*. This length of\n"
          "this array must match the size of those WSVs.\n"
          "\n"
          "Usage: Set by the user.\n"),
      GROUP("ArrayOfString")));

  wsv_data.push_back(WsvRecord(
      NAME("sideband_response"),
      DESCRIPTION(
          "Description of (mixer) sideband response.\n"
          "\n"
          "This variable describes the response of each sideband of a heterodyne\n"
          "receiver. The response is given as a GriddedField1, with frequency as the\n"
          "grid. The actual data describe the sideband filter function at each\n"
          "frequency grid point. An interpolation is applied to obtain the\n"
          "response for other frequencies.\n"
          "\n"
          "The frequency grid should be given in terms of IF, with end points\n"
          "symmetrically placed around zero. That is, the grid must contain\n"
          "both negative and positive values. The sideband response (after \n"
          "summation with *lo*) is not allowed to extend outside the range\n"
          "for which spectral data exist (normally determined by *f_grid*).\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Dimensions: \n"
          "\n"
          "-  GriddedField1:\n"
          "\n"
          "  -    Vector f_grid[N_f]\n"
          "  -    Vector data[N_f]\n"),
      GROUP("GriddedField1")));

  wsv_data.push_back(WsvRecord(
      NAME("sideband_response_multi"),
      DESCRIPTION(
          "Description of multiple (mixer) sideband responses.\n"
          "\n"
          "As *sideband_response* but describes an instrument with multiple\n"
          "mixers. An array element for each LO. The size of this variable and\n"
          "*lo_multi* shall match.\n"
          "\n"
          "Unit: Hz\n"
          "\n"
          "Usage: Set by the user.\n"),
      GROUP("ArrayOfGriddedField1")));

  wsv_data.push_back(WsvRecord(
      NAME("spectral_irradiance_field"),
      DESCRIPTION(
          "Spectral irradiance field.\n"
          "\n"
          "Spectral irradiance is the radiative power per unit area\n"
          "and unit frequency. The last dimension denotes the hemispheres.\n"
          "The first component denotes the downward direction and the second\n"
          "component denotes the upward direction.\n"
          "\n"
          "Units: W m^-2 Hz^-1\n"
          "\n"
          " Size: [ Nf, p_grid, lat_grid, lon_grid, 2 ]\n"),
      GROUP("Tensor5")));

  wsv_data.push_back(WsvRecord(
      NAME("spectral_radiance_field"),
      DESCRIPTION(
          "Spectral radiance field.\n"
          "\n"
          "This variable holds a calculation of the radiance field through\n"
          "the atmosphere, for the directions matching *za_grid* and *aa_grid*.\n"
          "\n"
          "Don't confuse this variable with *cloudbox_field*. That varinale also\n"
          "holds a field of spectral radiances, but is restricted to the cloud box.\n"
          "\n"
          "Units: W / (m^2 Hz sr)\n"
          "\n"
          " Size: [f_grid, p_grid,  lat_grid,  lon_grid,  za_grid, aa_grid, stokes_dim ]\n"
          "\n"
          "Note:\n"
          " For 1D, the size of the latitude, longitude and azimuth\n"
          " dimension (N_aa) are all 1.\n"),
      GROUP("Tensor7")));

  wsv_data.push_back(WsvRecord(
      NAME("specific_heat_capacity"),
      DESCRIPTION("Specific heat capacity.\n"
                  "\n"
                  "It is the heat capacity per unit \n"
                  "mass of a material.\n"
                  "\n"
                  "Units: K J^-1 kg^-1\n"
                  "\n"
                  "Size: [ p_grid, lat_grid,  lon_grid] \n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("specular_los"),
      DESCRIPTION(
          "The specular direction (for reflection by a flat surface).\n"
          "\n"
          "The specular direction as a standard line-of-sight vector, consisting\n"
          "of a zenith and azimuth angle (the later only for 3D).\n"
          "\n"
          "Units: degrees\n"
          "\n"
          "Size:  [ 1 or 2 ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("spt_calc_agenda"),
      DESCRIPTION(
          "Agenda calculating single scattering properties from the amplitude matrix.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(
      WsvRecord(NAME("stokes_dim"),
                DESCRIPTION("The dimensionality of the Stokes vector (1-4).\n"
                            "\n"
                            "Usage:      Set by the user.\n"),
                GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("suns_do"),
      DESCRIPTION("Flag to activate the sun(s).\n"),
      GROUP("Index"), Index{0}));

  wsv_data.push_back(WsvRecord(
      NAME("suns"),
      DESCRIPTION("Array of Sun.\n"
                  "\n"
                  "This variable describes a list of suns.\n"
                  "Each sun is described by a struct with its spectrum, radius,\n"
                  "distance from center of planet to center of sun,\n"
                  "temperature (if possible), latitude in the sky of the planet,\n"
                  "longitude in the sky of the planet and the type\n"),
      GROUP("ArrayOfSun"), ArrayOfSun{}));

  wsv_data.push_back(WsvRecord(
      NAME("stokes_rotation"),
      DESCRIPTION(
          "Rotation of the Stokes H and V directions.\n"
          "\n"
          "This variable allows to introduce a rotation of the Stokes coordinate\n"
          "system. Such a rotation could be needed to handle the scanning\n"
          "procedure of some instruments, such as AMSU-A. The variable is\n"
          "applied by the *sensor_responseStokesRotation* WSM.\n"
          "\n"
          "The rotation is given as an angle for each direction. In general, the\n"
          "number of rotations to be specified follows *sensor_response_dlos_grid*.\n"
          "In more detail, if no antenna is included or a 1D antenna is used, and\n"
          "the rotation is applied before the antenna is included in \n"
          "*sensor_response*, there should be one angle for each row of\n"
          "*mblock_dlos*. After inclusion of an antenna response, the relevant\n"
          "number of angles is determined by the rows of *antenna_dlos*.\n"
          "\n"
          "It is assumed that the rotation is common for all frequency elements.\n"
          "\n"
          "Units: degrees\n"
          "\n"
          "Size:  [ number of directions ]\n"
          "\n"
          "Usage: Set by the user.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("surface_complex_refr_index"),
      DESCRIPTION(
          "Complex refractive index of the surface, at a single point.\n"
          "\n"
          "See *complex_refr_index* for the expected format and how the data\n"
          "are treated.\n"),
      GROUP("GriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("surface_emission"),
      DESCRIPTION(
          "The emission from the surface.\n"
          "\n"
          "See specific methods generating *surface_emission* and the user\n"
          "guide for more information.\n"
          "\n"
          "Dimensions: [ f_grid, stokes_dim ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("surface_los"),
      DESCRIPTION(
          "Downwelling radiation directions to consider in surface reflection.\n"
          "\n"
          "The directions are given as a zenith and azimuth angle (the later\n"
          "only for 3D), following the definition of line-of-sights.\n"
          "\n"
          "Units: degrees\n"
          "\n"
          "Size:  [ any number, 1 or 2 ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("surface_normal"),
      DESCRIPTION(
          "The normal vector for a point at the surface.\n"
          "\n"
          "The vector is given as a zenith and azimuth (the later only for 3D)\n"
          "angle, following the definition of line-of-sights. For example,\n"
          "this vector is always [0] for 1D, as there is no surface topography\n"
          "for this atmospheric dimensionality.\n"
          "\n"
          "Units: degrees\n"
          "\n"
          "Size:  [ 1 or 2 ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("surface_props_data"),
      DESCRIPTION(
          "Various surface properties.\n"
          "\n"
          "A general container for passing data to surface methods. Each surface\n"
          "property shall be specified on the grid set by *lat_grid* and *lon_grid*.\n"
          "\n"
          "The properties are identified by the accompanying variable\n"
          "*surface_props_names*.\n"
          "\n"
          "Size:  [ number of props., lat_grid, lon_grid ]\n"),
      GROUP("Tensor3"), Tensor3{}));

  wsv_data.push_back(WsvRecord(
      NAME("surface_props_names"),
      DESCRIPTION(
          "Name on surface properties found in *surface_props_data*.\n"
          "\n"
          "Each string names a property in *surface_props_data*. The user is free\n"
          "to include data with any name, but the surface methods making use of\n"
          "*surface_props_data* expect data to be named in a specific way. See\n"
          "the documentation of each method for recognised choices.\n"
          "\n"
          "Size:  [ number of props. ]\n"),
      GROUP("ArrayOfString"), ArrayOfString{}));

  wsv_data.push_back(WsvRecord(
      NAME("surface_rmatrix"),
      DESCRIPTION(
          "The reflection coefficients for the directions given by\n"
          "*surface_los* to the direction of interest.\n"
          "\n"
          "The rows and columns of this tensor holds the reflection\n"
          "coefficient matrix for one frequency and one LOS. The reflection\n"
          "coefficients shall take into accound the angular weighting of the\n"
          "downwelling radiation.\n"
          "\n"
          "See specific methods generating *surface_rmatrix* and the user guide\n"
          "for more information.\n"
          "\n"
          "Usage:      Input to methods for *surface_rtprop_agenda*.\n"
          "\n"
          "Units:      -\n"
          "\n"
          "Dimensions: [ surface_los, f_grid, stokes_dim, stokes_dim ]\n"),
      GROUP("Tensor4")));

  wsv_data.push_back(WsvRecord(
      NAME("surface_rtprop_agenda"),
      DESCRIPTION("Agenda providing radiative properties of the surface.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("surface_rtprop_agenda_array"),
      DESCRIPTION(
          "Description of surface radiative properties, for each surface type.\n"),
      GROUP("ArrayOfAgenda")));

  wsv_data.push_back(WsvRecord(
      NAME("surface_skin_t"),
      DESCRIPTION(
          "Surface skin temperature.\n"
          "\n"
          "This temperature shall be selected considering the radiative\n"
          "properties of the surface, and can differ from the \"bulk\"\n"
          "temperature.\n"
          "\n"
          "Usage:   Input to methods for *surface_rtprop_agenda*.\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("surface_reflectivity"),
      DESCRIPTION(
          "Surface reflectivity, for a given position and angle.\n"
          "\n"
          "This variable describes the surface reflectivity at one position\n"
          "and one incidence angle. It works as *surface_scalar_reflectivity*\n"
          "but is also defined for vector radiative transfer.\n"
          "\n"
          "The first dimension of the variable shall either match *f_grid* or\n"
          "be 1. The later case is interpreted as the reflectivity is the same\n"
          "for all frequencies.\n"
          "\n"
          "Usage:   Input to some surface properties methods.\n"
          "\n"
          "Dimensions: [ f_grid or 1, stokes_dim, stokes_dim]\n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("surface_rv_rh"),
      DESCRIPTION(
          "Surface reflectivity, described by rv and rh (power) reflectivities.\n"
          "\n"
          "This variable describes the surface reflectivity at one position\n"
          "and one incidence angle. For this position and angle, one or multiple\n"
          "combinations of rv and rh are specified, where rv and rh are the\n"
          "reflectivity for vertical and horizontal polarisation, respectively.\n"
          "\n"
          "This matrix shall always have two columns, where the first column\n"
          "holds rv values, and the second column rh. It is up to the user to\n"
          "make sure that data are put into the correct column, this can not\n"
          "be checked bu the methods using this WSV.\n"
          "\n"
          "The number of rows shall either match *f_grid* or be 1. The later case\n"
          "is interpreted as the reflectivities are the same for all frequencies.\n"
          "\n"
          "Usage:   Input to some surface properties methods.\n"
          "\n"
          "Dimensions: [ f_grid or 1, 2]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("surface_scalar_reflectivity"),
      DESCRIPTION(
          "Surface reflectivity, assuming it can be described as a scalar value.\n"
          "\n"
          "This variable describes the surface reflectivity at one position\n"
          "and one incidence angle. For this position and angle, one or multiple\n"
          "scalar reflectivities are specified.\n"
          "\n"
          "The length of the vector shall either match *f_grid* or be 1. The \n"
          "later case is interpreted as the reflectivity is the same for all\n"
          "frequencies (ie. matches a constant vector).\n"
          "\n"
          "Usage:   Input to some surface properties methods.\n"
          "\n"
          "Dimensions: [ f_grid or 1]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("surface_type_mask"),
      DESCRIPTION(
          "Classification of the surface using a type coding.\n"
          "\n"
          "There is no fixed type coding, it is up to the user to set up\n"
          "a system consistent with *surface_rtprop_agenda_array*. A value\n"
          "of 0 in *surface_type_mask* means that element 0 in the agenda\n"
          "array is valid for that position etc.\n"
          "\n"
          "Dimensions: \n"
          "\n"
          "-  GriddedField2:\n"
          "\n"
          "  -    Vector Latitude [N_lat]\n"
          "  -    Vector Longitude [N_lon]\n"
          "  -    Matrix data [N_lat][N_lon]\n"),
      GROUP("GriddedField2")));

  wsv_data.push_back(
      WsvRecord(NAME("surface_type_mix"),
          DESCRIPTION(
            "Gives the fraction of different surface types.\n"
            "\n"
            "For cases when the surface RT properties are taken from\n"
            "*surface_rtprop_agenda_array*, this variable specifies to\n"
            "what extent each surface type has contributed to the surface\n"
            "RT variables, such as *surface_emission* and *surface_skin_t*.\n" 
            "\n"
            "The length of this vector follows *surface_rtprop_agenda_array*\n"
            "and the sum of the elements is 1. The first element in the\n" 
            "vector matches the first agenda element, and so on."),
          GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("telsem_atlases"),
      DESCRIPTION(
          "TELSEM 2 emissivity atlases.\n"
          "\n"
          "Array should be filled with 12\n"
          "atlases, one for each month. Index 0 is January, index 11 December.\n"
          ""),
      GROUP("ArrayOfTelsemAtlas")));

  wsv_data.push_back(WsvRecord(
      NAME("tessem_neth"),
      DESCRIPTION(
          //FIXMEDOC Add more documentation?
          "TESSEM2 neural network parameters for horizontal polarization.\n"),
      GROUP("TessemNN")));

  wsv_data.push_back(WsvRecord(
      NAME("tessem_netv"),
      DESCRIPTION(
          //FIXMEDOC Add more documentation?
          "TESSEM2 neural network parameters for vertical polarization.\n"),
      GROUP("TessemNN")));

  wsv_data.push_back(
      WsvRecord(NAME("test_agenda"),
                DESCRIPTION("A dummy agenda for testing purposes.\n"
                            "\n"
                            "Only used for testing by developers.\n"),
                GROUP("Agenda")));

  wsv_data.push_back(
      WsvRecord(NAME("test_agenda_array"),
                DESCRIPTION("Array of agenda for TestArrayOfAgenda case.\n"
                            "\n"
                            "Only used for testing by developers.\n"),
                GROUP("ArrayOfAgenda")));

  wsv_data.push_back(WsvRecord(
      NAME("time"),
      DESCRIPTION("A UTC time point.\n"),
      GROUP("Time")));

  wsv_data.push_back(WsvRecord(
      NAME("timer"),
      DESCRIPTION("Stores the starting time for time measurements.\n"),
      GROUP("Timer")));

  wsv_data.push_back(WsvRecord(
      NAME("time_grid"),
      DESCRIPTION("A grid of times.  Should be increasing\n"),
      GROUP("ArrayOfTime")));

  wsv_data.push_back(WsvRecord(
      NAME("time_stamps"),
      DESCRIPTION("A set of times.  Can be in random order\n"),
      GROUP("ArrayOfTime")));

  wsv_data.push_back(WsvRecord(
      NAME("transmitter_pos"),
      DESCRIPTION(
          "Transmitter positions.\n"
          "\n"
          "Used for radio link calculations and gives then the position of the\n"
          "transmitting device. The corresponding positions of the receiver are\n"
          "given by *sensor_pos*. The number of rows in *transmitter_pos* and\n"
          "*sensor_pos* must be equal.\n"
          "\n"
          "This WSV is also defined as *sensor_pos* regarding the content of the\n"
          "columns, accepted range for latitudes etc. With one exception, this\n"
          "WSV is demanded to have two columns also for 1D. The additional\n"
          "second value is the angular distance between the transmitter and the\n"
          "reciver. This angle is defined as \"latitude\" for 2D, with the\n"
          "sensor fixed at the angle of 0 degree.\n"
          "\n"
          "Each row this matrix defines *rte_pos2* for the measurement block,\n"
          "exactly as *sensor_pos* is translated to *rte_pos*.\n"
          "\n"
          "If no transmitter is involved in the calculations, the variable can\n"
          "be set to be empty.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  [ m, degrees, degrees ]\n"),
      GROUP("Matrix"), Matrix{}));

  wsv_data.push_back(WsvRecord(
      NAME("t_field"),
      DESCRIPTION(
          "The field of atmospheric temperatures.\n"
          "\n"
          "This variable gives the atmospheric temperature at each crossing of\n"
          "the pressure, latitude and longitude grids.\n"
          "\n"
          "The temperature for a point between the grid crossings is obtained \n"
          "by (multi-)linear interpolation of the *t_field*.\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage:      Output of *AtmFieldsCalc*.\n"
          "\n"
          "Unit:       K\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ]\n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("nlte_field"),
      DESCRIPTION(
          "The field of NLTE temperatures and/or ratios.\n"
          "\n"
          "This variable gives the NLTE temperature/ratio at each crossing of\n"
          "the pressure, latitude and longitude grids.  The size of the\n"
          "array is the number of NLTE levels in all molecules.\n"
          "\n"
          "The temperature/ratio for a point between the grid crossings is obtained \n"
          "by (multi-)linear interpolation of the *nlte_field*.\n"
          "\n"
          "There are two types of NLTE computations available in ARTS.  One from\n"
          "giving excitiation temperatures that makes the absorption/emission diverge\n"
          "from LTE.  The other is to use the absolute ratios of upper-to-lower states at\n"
          "the levels of interest.\n"
          ""
          "\n"
          "Units:       [ K or \% ]]\n"
          "\n"
          "Dimensions: [ NLTE levels, p_grid, lat_grid, lon_grid ] or [ 0, 0, 0, 0 ]\n"),
      GROUP("EnergyLevelMap"), EnergyLevelMap{}));

  wsv_data.push_back(WsvRecord(
      NAME("t_field_raw"),
      DESCRIPTION(
          "Raw data for atmospheric temperatures.\n"
          "\n"
          "This variable gives the atmospheric temperature as stored in the \n"
          "database for the atmospheric scenarios.\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage: Set by the user by choosing a climatology.\n"
          "\n"
          "Unit:  K\n"
          "\n"
          "Size:\n"
          "\n"
          "-  GriddedField3\n"
          "\n"
          "  -     [N_p]\n"
          "  -     [N_lat]\n"
          "  -     [N_lon]\n"
          "  -     [N_p, N_lat, N_lon] \n"),
      GROUP("GriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("nlte_field_raw"),
      ("Raw data for NLTE temperatures and/or ratios.\n"
       "\n"
       "This variable gives the NLTE temperature/ratio as stored in the \n"
       "database for the atmospheric scenarios.\n"
       "\n"
       "See further the ARTS user guide (AUG). Use the index to find where\n"
       "this variable is discussed. The variable is listed as a subentry to\n"
       "\"workspace variables\".\n"
       "\n"
       "Usage: Set by the user by choosing a climatology.\n"
       "\n"
       "Unit:  K\n"
       "\n"
       "Size:\n"
       "\n"
       "- ArrayOfGriddedField3 \n"
       "\n"
       "  -     [NLTE levels] or [ 0 ]\n"
       "\n"
       "    -    [N_p] \n"
       "    -    [N_lat] \n"
       "    -    [N_lon] \n"
       "    -    [N_p, N_lat, N_lon] \n"),
      GROUP("ArrayOfGriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("t_surface"),
      DESCRIPTION(
          "The surface temperature.\n"
          "\n"
          "This variable holds the temperature of the surface at each latitude\n"
          "and longitude grid crossing. The normal case should be that this \n"
          "temperature field is interpolated to obtain *surface_skin_t*.\n"
          "Accordingly, for 1D cases it could be a better idea to specify\n"
          "*surface_skin_t* directly.\n"
          "\n"
          "These temperature shall be selected considering the radiative\n"
          "properties of the surface, and can differ from the \"bulk\"\n"
          "temperatures.\n"
          "\n"
          "Usage:      Set by user.\n"
          "\n"
          "Unit:       K\n"
          "\n"
          "Dimensions: [ lat_grid, lon_grid ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("verbosity"),
      DESCRIPTION(
          "ARTS verbosity.\n"
          "\n"
          "The verbosity variable is implicitly passed to all workspace methods.\n"
          "It can be used to dynamically control the reporting level during\n"
          "runtime.\n"
          "\n"
          "Usage:    Set by user.\n"
          "\n"
          "See also:\n"
          "          *verbosityInit*\n"
          "          *verbositySet*\n"
          "          *verbositySetAgenda*\n"
          "          *verbositySetScreen*\n"
          "          *verbositySetFile*\n"),
      GROUP("Verbosity")));

  wsv_data.push_back(WsvRecord(
      NAME("vmr_field"),
      DESCRIPTION(
          "VMR field.\n"
          "\n"
          "This variable gives the volume mixing ratio of the chosen gaseous \n"
          "species as a function of p_grid, lat_grid, lon_grid. \n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Units: [ Differ between the elements, can be VMR, kg/m3 or #/m3. ]\n"
          "\n"
          "Dimensions: [species, p_grid, lat_grid, lon_grid]\n"),
      GROUP("Tensor4")));

  wsv_data.push_back(WsvRecord(
      NAME("vmr_field_raw"),
      DESCRIPTION(
          "VMR data for the chosen gaseous species.\n"
          "\n"
          "This variable contains the volume mixing ratios (VMR) for all \n"
          "chosen gaseous species. It includes the grids corresponding to the \n"
          "grids in the database. \n"
          "*vmr_field_raw* is an Array of Array of Tensor3. It contains one \n"
          "gridded field for each species which contains the data and \n"
          "also the grids.\n"
          "For the calculation the data is \n"
          "interpolated on *p_grid*, *lat_grid* and *lon_grid*\n"
          "\n"
          "Usage:\n"
          "\n"
          "- Output of *AtmRawRead*.\n"
          "- Input to *AtmFieldsCalc*.\n"
          "\n"
          "Unit:  absolute number\n"
          "\n"
          "Size:\n"
          "\n"
          "- Array[number of absorption species]\n"
          "\n"
          "  -     GriddedField3\n"
          "\n"
          "    -    [N_p] \n"
          "    -    [N_lat] \n"
          "    -    [N_lon] \n"
          "    -    [N_p, N_lat, N_lon] \n"),
      GROUP("ArrayOfGriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("water_p_eq_agenda"),
      DESCRIPTION("Agenda to calculate the saturation pressure of water.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("water_p_eq_field"),
      DESCRIPTION(
          "The field of water saturation pressure.\n"
          "\n"
          "This variable holds the saturation pressure of water at each crossing of\n"
          "the pressure, latitude and longitude grids.\n"
          "\n"
          "Unit:       Pa\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ]\n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("wigner_initialized"),
      DESCRIPTION(
          "Indicates if the wigner tables are initialized.\n"
          "If they are not, computations will be aborted.\n"
          "\n"
          "Will hold the value of provided maximum factorial value\n"
          "\n"
          "The developer should always test this variable in functions\n"
          "that might require computing wigner symbols because the error\n"
          "handling is otherwise offloaded to third party software...\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("wind_u_field"),
      DESCRIPTION(
          "Zonal component of the wind field.\n"
          "\n"
          "The East-West wind component. Air moving towards higher\n"
          "longitudes is a positive wind. This wind causes no Doppler shift\n"
          "for 1D and 2D simulations.\n"
          "\n"
          "Can be set to be empty, which is interpreted as zero wind speed\n"
          "everywhere.\n"
          "\n"
          "Unit:       m/s\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ]  or [ 0 0 0 ].\n"),
      GROUP("Tensor3"), Tensor3{}));

  wsv_data.push_back(WsvRecord(
      NAME("wind_u_field_raw"),
      DESCRIPTION(
          "Raw zonal component of the wind field.\n"
          "\n"
          "The East-West wind component. Air moving towards higher\n"
          "longitudes is a positive wind. This wind causes no Doppler shift\n"
          "for 1D and 2D simulations.\n"
          "\n"
          "Can be set to be empty, which is interpreted as zero wind speed\n"
          "everywhere.\n"
          "\n"
          "Unit:       m/s\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ].\n"),
      GROUP("GriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("wind_v_field"),
      DESCRIPTION(
          "Meridional component of the magnetic field.\n"
          "\n"
          "The North-South wind component. Air moving towards higher\n"
          "latitudes is a positive wind.\n"
          "\n"
          "Can be set to be empty, which is interpreted as zero wind speed\n"
          "everywhere.\n"
          "\n"
          "Unit:       m/s\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ] or [ 0 0 0 ]\n"),
      GROUP("Tensor3"), Tensor3{}));

  wsv_data.push_back(WsvRecord(
      NAME("wind_v_field_raw"),
      DESCRIPTION(
          "Raw meridional component of the magnetic field.\n"
          "\n"
          "The North-South wind component. Air moving towards higher\n"
          "latitudes is a positive wind.\n"
          "\n"
          "Can be set to be empty, which is interpreted as zero wind speed\n"
          "everywhere.\n"
          "\n"
          "Unit:       m/s\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ]\n"),
      GROUP("GriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("wind_w_field"),
      DESCRIPTION(
          "Vertical wind component field.\n"
          "\n"
          "Upward moving air corresponds to a positive wind speed.\n"
          "\n"
          "Can be set to be empty, which is interpreted as zero wind speed\n"
          "everywhere.\n"
          "\n"
          "Unit:       m/s\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ] or [ 0 0 0 ]\n"),
      GROUP("Tensor3"), Tensor3{}));

  wsv_data.push_back(WsvRecord(
      NAME("wind_w_field_raw"),
      DESCRIPTION(
          "Raw vertical wind component field.\n"
          "\n"
          "Upward moving air corresponds to a positive wind speed.\n"
          "\n"
          "Can be set to be empty, which is interpreted as zero wind speed\n"
          "everywhere.\n"
          "\n"
          "Unit:       m/s\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ]\n"),
      GROUP("GriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("wmrf_channels"),
      DESCRIPTION(
          "Channel selection for WMRF fast calculation.\n"
          "\n"
          "This variable can be used to select one or several instrument channels\n"
          "from the list of all possible channels. Zero-based indexing is used, so\n"
          "Channel 0 is the first instrument channel!\n"),
      GROUP("ArrayOfIndex")));

  wsv_data.push_back(WsvRecord(
      NAME("wmrf_weights"),
      DESCRIPTION(
          "The weights for a WMRF fast calculation.\n"
          "\n"
          "Weights are stored in a sparse matrix. This can be used as a\n"
          "sensor_response matrix.\n"
          "\n"
          "The dimension of the matrix is (nchan, nfreq), where nchan\n"
          "is the number of instrument channels and nfreq is the number\n"
          "of monochromatic frequencies.\n"),
      GROUP("Sparse")));

  wsv_data.push_back(WsvRecord(
      NAME("xml_output_type"),
      DESCRIPTION(
          "Flag to determine whether XML output shall be binary or ascii.\n"
          "\n"
          "This flag has to be set using the workspace method\n"
          "*output_file_formatSetAscii* or *output_file_formatSetBinary*.\n"
          "One of these methods MUST be called before writing the first\n"
          "output file.\n"
          "\n"
          "Usage: Set by user.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("x"),
      DESCRIPTION(
          "The state vector.\n"
          "\n"
          "This WSV matches directly the x-vector in the formalism by C.D. Rodgers.\n"
          "\n"
          "Inside *x*, the elements matching one retrieval quantity, such as\n"
          "atmospheric temperatures, are kept together. That is, each retrieval\n"
          "quantity covers a continuous range inside *x*. The start and index of\n"
          "these ranges can be deduced by *jacobian_quantities* (see function(s)\n"
          "inside jacobian.cc for details).\n"
          "\n"
          "The order of elements inside each retrieval quantity should be clarified\n"
          "by corresponding \"adding\" method, i.e. *jacobianAddTemperature* for\n"
          "atmospheric temperatures. The general rule is that data are sorted from\n"
          "left to right with respect to the order in the corresponding WSV. For\n"
          "example, inside *x* atmospheric data are stored with pressure as inner-\n"
          "most loop, followed by latitude and longitude as outermost loop.\n"
          "\n"
          "Usage: Used by inversion methods.\n"
          "\n"
          "Unit:  Varies, follows unit of selected retrieval quantities.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("xa"),
      DESCRIPTION(
          "The a priori state vector.\n"
          "\n"
          "This WSV matches directly the x_a-vector in the formalism by C.D. Rodgers.\n"
          "\n"
          "Usage: Used by inversion methods.\n"
          "\n"
          "Unit:  Varies, follows unit of selected retrieval quantities.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("y"),
      DESCRIPTION(
          "The measurement vector.\n"
          "\n"
          "This vector holds radiances averaged in frequency and spatially,\n"
          "and can contain many spectra appended. That is, this WSV matches\n"
          "directly the y-vector in the formalism by C.D. Rodgers.\n"
          "\n"
          "The polarisation, frequency, position and line-of-sight associated\n"
          "with each element in *y* are given by *y_pol*, *y_f*, *y_pos* and\n"
          "*y_los*. For monochromatic pencil beam radiances, data are sorted\n"
          "in the following way, from the innermost to the outermost loop\n"
          "\n"
          "-   Stokes\n"
          "-   Frequency\n"
          "-   LOS inside the measurement block\n"
          "-   Measurement block\n"
          "\n"
          "With sensor response included, the order can be differ. As output\n"
          "of *yRadar*, the order will also be different.\n"
          "\n"
          "Usage: Output from radiative transfer calculations considering sensor response.\n"
          "\n"
          "Unit:  Undefined. Possibilities include: K, W/(m^2 Hz sr) and optical thickness.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("y_aux"),
      DESCRIPTION(
          "Data auxilary to *y*.\n"
          "\n"
          "Different data beside the direct result of the radiative transfer\n"
          "calculations can be obtained by this variable. These auxilary data\n"
          "are selected by *iy_aux_vars*.\n"
          "\n"
          "In contrast to *iy_aux*, this variable can only hold quantities such\n"
          "as optical depth, and other quantites that could be the result\n"
          "of a complete radiative transfer calculation. The data are weighted\n"
          "with sensor properties in the same way as for *y*.\n"
          "\n"
          "See also *iy_aux_vars*.\n"
          "\n"
          "Usage:      Output of *yCalc*.\n"
          "\n"
          "Dimensions: [quantity][ element of y ]\n"),
      GROUP("ArrayOfVector")));

  wsv_data.push_back(WsvRecord(
      NAME("y_baseline"),
      DESCRIPTION(
          "The baseline of *y*.\n"
          "\n"
          "In retrieval \"lingo\", the baseline is an addiative disturbance of\n"
          "the measured spectrum. That is, it can be seen as a shift (from zero)\n"
          "of measurement. Reflections inside microwave receivers is one source to\n"
          "a baseline off-set.\n"
          "\n"
          "So far there is no module in ARTS that actually tries to physically model\n"
          "any baseline effect. *y_baseline* is just used as a pure fitting parameter\n"
          "in retrievals. One example on method to include a baseline fit is \n"
          "*jacobianAddPolyfit*.\n"
          "\n"
          "If the baseline is totally constant, it is allowed to set *y_baseline*\n"
          "to have length one, with this element set to the baseline value.\n"
          "\n"
          "Usage: Output of retrievals.\n"
          "\n"
          "Unit:  Same as applied for *y*.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("y_f"),
      DESCRIPTION(
          "The frequencies associated with *y*.\n"
          "\n"
          "A value is returned for each element of *y*. Depending on the sensor\n"
          "set-up and number of measurement blocks, this can be a copy of\n"
          "*sensor_response_f*, sveral copies of this vector appended, or some\n"
          "other frequenices.\n"
          "\n"
          "Don't confuse this variable with *yf*.\n"
          "\n"
          "Usage: Output from radiative transfer calculations considering sensor response.\n"
          "\n"
          "Unit:  [ Hz ]\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("y_geo"),
      DESCRIPTION(
          "The geo-position assigned to each element of  *y*.\n"
          "\n"
          "The columns of this matrix matches the elements of *geo_pos*.\n"
          "\n"
          "Unit:  [ m, deg, deg, deg, deg ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("y_geo_series"),
      DESCRIPTION(
          "The geo-positioning assigned to each row of *y_series*.\n"
          "\n"
          "All channels are assumed to have the same geo-position.\n"
          "\n"
          "Otherwise as *y_geo*.\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("y_geo_swath"),
      DESCRIPTION(
          "The geo-positioning assigned to each pixel of *y_swath*.\n"
          "\n"
          "All channels are assumed to have the same geo-position.\n"
          "\n"
          "Otherwise as *y_geo*.\n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("y_los"),
      DESCRIPTION(
          "The line-of-sights associated with *y*.\n"
          "\n"
          "Definition of angles matches *sensor_los* (such as first column holds\n"
          "zenith angles), but gives actual observed LOS. That is, the values of\n"
          "both *sensor_los* and *antenna_dlos* are considered. Data are provided\n"
          "for each element of *y*, following y_f, and the number of rows equals\n"
          "the length of *y*.\n"
          "\n"
          "Usage: Output from radiative transfer calculations considering sensor response.\n"
          "\n"
          "Unit:  [ degrees, degrees ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("y_pol"),
      DESCRIPTION(
          "The polarisation states associated with *y*.\n"
          "\n"
          "Data are provided for each element of *y*, following y_f, and the\n"
          "length of this variable and *y* is equal.\n"
          "\n"
          "See *instrument_pol* for coding of polarisation components.\n"
          "\n"
          "Usage: Output from radiative transfer calculations considering sensor response.\n"
          "\n"
          "Unit:  [ - ]\n"),
      GROUP("ArrayOfIndex")));

  wsv_data.push_back(WsvRecord(
      NAME("y_pos"),
      DESCRIPTION(
          "The sensor positions associated with *y*.\n"
          "\n"
          "Definition of positions matches *sensor_pos* (such as first column\n"
          "holds the altitude). Data are provided for each element of *y*,\n"
          "following y_f, and the number of rows equals the length of *y*.\n"
          "\n"
          "Usage: Output from radiative transfer calculations considering sensor response.\n"
          "\n"
          "Unit:  [ m, deg, deg ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("y_series"),
      DESCRIPTION(
          "Two-dimensional version of the measurement vector.\n"
          "\n"
          "This WSV can be used for storing *y* reshaped when all measurement\n"
          "blocks have the same set of channels.\n"
          "\n"
          "Dimesion:  [ position, channel ]\n"),
      GROUP("Matrix")));

  wsv_data.push_back(WsvRecord(
      NAME("y_swath"),
      DESCRIPTION(
          "Three-dimensional version of the measurement vector.\n"
          "\n"
          "This WSV can be used for storing *y* reshaped when all measurement\n"
          "blocks have the same set of channels, and that the data constitutes\n"
          "a part of a swath.\n"
          "\n"
          "Dimesion:  [ scan, pixel, channel ]\n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("yb"),
      DESCRIPTION(
          "The measurement vector for a single measurement block.\n"
          "\n"
          "Exactly as *y*, but holds data only for a single measurement block.\n"
          "\n"
          "Usage: Used internally.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("ybatch"),
      DESCRIPTION(
          "Batch of spectra.\n"
          "\n"
          "Each element of *ybatch* corresponds to a spectrum vector *y*. \n"
          "See further *ybatchCalc*.\n"
          "\n"
          "Usage: Most commonly produced by *ybatchCalc*.\n"
          "\n"
          "Unit:  Undefined. Possibilities include: K, W/(m^2 Hz sr) and optical thickness.\n"
          "\n"
          "Dimensions: Number of array elements equals number of batch cases,\n"
          "Vectors have length(y)\n"),
      GROUP("ArrayOfVector")));

  wsv_data.push_back(WsvRecord(
      NAME("ybatch_aux"),
      DESCRIPTION(
          "Data auxilary to *ybatch*.\n"
          "\n"
          "Each element of *ybatch_aux* corresponds to a auxiliary data *y_aux*. \n"
          "See further *y_aux* and *ybatchCalc*.\n"
          "\n"
          "Usage: Most commonly produced by *ybatchCalc*.\n"
          "\n"
          "Dimensions: Number of array elements equals number of batch cases,\n"),
      GROUP("ArrayOfArrayOfVector")));

  wsv_data.push_back(WsvRecord(
      NAME("ybatch_calc_agenda"),
      DESCRIPTION(
          "Agenda defining the calculations to perform for each batch case.\n"),
      GROUP("Agenda")));

  wsv_data.push_back(WsvRecord(
      NAME("ybatch_index"),
      DESCRIPTION("Index of batch case.\n"
                  "\n"
                  "See further *ybatchCalc*.\n"
                  "\n"
                  "Usage: Set by *ybatchCalc*, for communication with\n"
                  "*ybatch_calc_agenda*.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("ybatch_corr"),
      DESCRIPTION(
          "Correction terms for *ybatch*.\n"
          "\n"
          "Dimensions: Number of array elements equals number of batch cases,\n"
          "Vectors have length depending on correction method\n"),
      GROUP("ArrayOfVector")));

  wsv_data.push_back(WsvRecord(
      NAME("ybatch_jacobians"),
      DESCRIPTION("All the Jacobians associated with ybatch.\n"
                  "\n"
                  "The batch index here is the array dimension.\n"
                  "\n"
                  "Usage: Most commonly produced by *ybatch*.\n"
                  "\n"
                  "Unit:  Depends on unit of y and on Jacobian type.\n"
                  "\n"
                  "Dimensions:\n"
                  "\n"
                  "- [number of batch cases]\n"
                  "\n"
                  "  -           (length(y),\n"
                  "  -           number of retrieval quantities and grids)\n"),
      GROUP("ArrayOfMatrix")));

  wsv_data.push_back(
      WsvRecord(NAME("ybatch_n"),
                DESCRIPTION("Number of batch cases for *ybatchCalc*.\n"
                            "\n"
                            "See further *ybatchCalc*.\n"
                            "\n"
                            "Usage: Input to *ybatchCalc*.\n"),
                GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("ybatch_start"),
      DESCRIPTION("Start index for *ybatchCalc*.\n"
                  "\n"
                  "This is set to a default of zero.\n"
                  "\n"
                  "See further *ybatchCalc*.\n"
                  "\n"
                  "Usage: Input to *ybatchCalc*.\n"),
      GROUP("Index"), Index{0}));

  wsv_data.push_back(WsvRecord(
      NAME("yf"),
      DESCRIPTION(
          "A fitted measurement vector.\n"
          "\n"
          "This WSV is the measurement vector matching the retrieved state, i.e.\n"
          "the spectrum of the fit.\n"
          "\n"
          "Don't confuse this variable with *y_f*.\n"
          "\n"
          "Usage: Output from inversion methods.\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("za_grid"),
      DESCRIPTION(
          "Zenith angle grid.\n"
          "\n"
          "The zenith angle grid, on which the *cloudbox_field* is stored. \n"
          "This grid is used for RT calculations inside the cloudbox, therefore\n"
          "the grid has to be defined if the cloudbox is activated by the flag\n"
          "*cloudbox_on*. Furthermore the zenith angle grid is also used for RT\n"
          "calculations of clear-sky *spectral_radiance_field*. \n"
          "The grid must be sorted in increasing order, with no repetitions.\n"
          "\n"
          "Usage:      Set by the user.\n"
          "\n"
          "Unit:       degrees \n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("za_grid_weights"),
      DESCRIPTION("Zenith angle integration weights.\n"
          "\n"
          "The integration weight are needed for calculation of radiation fluxes\n"
          "\n"
          "Unit:  unitless\n"),
      GROUP("Vector")));

  wsv_data.push_back(WsvRecord(
      NAME("za_index"),
      DESCRIPTION(
          "Zenith angle index for scattering calculations.\n"
          " \n"
          "This variable is used internally in WSMs for computing scattering \n"
          "properties. \n"
          "\n"
          "Usage:    Input to the agendas *spt_calc_agenda*, *pha_mat_spt_agenda*.\n"),
      GROUP("Index")));

  wsv_data.push_back(WsvRecord(
      NAME("z_field"),
      DESCRIPTION(
          "The field of geometrical altitudes.\n"
          "\n"
          "This variable gives the geometrical altitude, above the ellipsoid, of\n"
          "each crossing of the pressure, latitude and longitude grids. For 1D\n"
          "cases the altitudes give the geometrical position of the pressure\n"
          "levels.\n"
          "\n"
          "For each geographical position (lat,lon), the values must be sorted\n"
          "in increasing order, with no repetitions. Otherwise the altitudes\n"
          "can be set to arbitrary values. Hydrostatic equilibrium is not\n"
          "applied automatically. If hydrostatic equilibrium applies, *z_field*\n"
          "must be set by a method ensuring that this criterium is fulfilled.\n"
          "\n"
          "The radius (from the coordinate centre) for a point between the grid\n"
          "crossings is obtained by a (multi-)linear interpolation of the sum\n"
          "of the ellipsoid radius and *z_field*.\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage:      Output of *AtmFieldsCalc*\n"
          "\n"
          "Unit:       m\n"
          "\n"
          "Dimensions: [ p_grid, lat_grid, lon_grid ]\n"),
      GROUP("Tensor3")));

  wsv_data.push_back(WsvRecord(
      NAME("z_field_raw"),
      DESCRIPTION(
          "Raw data for geometrical altitudes.\n"
          "\n"
          "This variable gives the geometrical altitudes as stored in the \n"
          "database for atmospheric scenarios.\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage: Set by the user by choosing a climatology.\n"
          "\n"
          "Unit:  K\n"
          "\n"
          "Size:\n"
          "\n"
          "-  GriddedField3\n"
          "\n"
          "  -     [N_p]\n"
          "  -     [N_lat]\n"
          "  -     [N_lon]\n"
          "  -     [N_p, N_lat, N_lon]\n"),
      GROUP("GriddedField3")));

  wsv_data.push_back(WsvRecord(
      NAME("z_hse_accuracy"),
      DESCRIPTION(
          "Minimum accuracy for calculation of hydrostatic equilibrium.\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  m\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("z_sensor"),
      DESCRIPTION(
          "The altitude of the sensor.\n"
          "\n"
          "Please note that the sensor altitude actaully applied is in general\n"
          "specified by *sensor_pos*. This WSV is only a help, to set other\n"
          "workspace variables and to call methods in a consistent manner\n"
          "\n"
          "Usage: Set by the user.\n"
          "\n"
          "Unit:  m\n"),
      GROUP("Numeric")));

  wsv_data.push_back(WsvRecord(
      NAME("z_surface"),
      DESCRIPTION(
          "The surface altitude.\n"
          "\n"
          "This variable defines the shape of the surface, by giving the\n"
          "geometrical altitude above the geiod for each crossing of the \n"
          "latitude and longitude grids. Any shape of the surface is accepted.\n"
          "No gap between the surface and the lowermost pressure level is \n"
          "allowed.\n"
          "\n"
          "The radius (from the coordinate centre) for a point between the grid\n"
          "crossings is obtained by a linear (1D) or bi-linear (2D) \n"
          "interpolation of the sum of the ellipsoid radius and *z_surface*.\n"
          "That is, the radius for the surface is assumed to vary linear along \n"
          "the latitudes and longitudes in *lat_grid* and *lon_grid*.\n"
          "\n"
          "See further the ARTS user guide (AUG). Use the index to find where\n"
          "this variable is discussed. The variable is listed as a subentry to\n"
          "\"workspace variables\".\n"
          "\n"
          "Usage:      Set by user.\n"
          "\n"
          "Unit:       m\n"
          "\n"
          "Dimensions: [ lat_grid, lon_grid ]\n"),
      GROUP("Matrix")));

  std::sort(wsv_data.begin(), wsv_data.end(), [](auto& a, auto& b) {
    return a.Name() < b.Name();
  });
}

void define_wsv_map() {
  for (Index i = 0; i < global_data::wsv_data.nelem(); ++i) {
    global_data::WsvMap[global_data::wsv_data[i].Name()] = i;
  }
}
