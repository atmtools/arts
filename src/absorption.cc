// Physical absorption routines. The absorption workspace methods are
// in file m_abs.cc 

#include "arts.h"
#include "make_array.h"
#include "absorption.h"

// Some #defines and typedefs to make the records better readable:
#define NAME(x) x
#define DEGFR(x) x
#define ISOTOPES     make_array<IsotopeRecord>
#define REC	     IsotopeRecord
#define TAGS	     make_array<int>


/** The lookup information for all the different species. */
ARRAY<SpeciesRecord> species_data;

/** The map associated with species_data. */
std::map<string, size_t> SpeciesMap;

void define_species_data()
{
  // Initialize to zero, just in case:
  species_data.clear();

  /* Here's an empty template record entry:

  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O"),
	DEGFR(3),
	ISOTOPES
	(//   Name,	Abundance,	Mass,	MY-tag, HY-tag, JPL-tag
	 //		|		|	|	|	|
	 REC( ""	,		,	,	,	,TAGS() ),
	 ) ) );

  For good looks, keep the commas on the marks!

  */


  // H2O
  // Source: 
  // Isotope ratios: HITRAN92 tables.3
  //		     ISOTOPERATIO for double isotopic species (not in
  //		     HITRAN) are products of individual ratios.	 
  species_data.push_back
    ( SpeciesRecord
      ( NAME("H2O"),
	DEGFR(3),
	ISOTOPES
	(//   Name,	Abundance,	Mass,	MY-tag, HY-tag, JPL-tag
	 //		|		|	|	|	|
	 REC( "161"	,0.997317	,18.	,11	,11	,TAGS(18003, 18005) ),
	 REC( "181"	,0.00199983	,20.	,12	,12	,TAGS(20003) ),
	 REC( "171"	,0.000372	,19.	,13	,13	,TAGS(19003) ),
	 REC( "162"	,0.00031069	,19.	,14	,14	,TAGS(19002) ),
	 REC( "182"	,621.327e-9	,21.	,-1	,-1	,TAGS(21001) ),
	 REC( "262"	,96.528e-9	,20.	,-1	,-1	,TAGS(20001) )
	 ) ) );

  // CO2 
  // (missing in JPL)
  // Source: 
  // Degrees of freedom from Schanda:`Physical Fundamentals of Remote Sensing'
  // FIXME: Did I really enter the degees of freedom? Check again in Schanda (Axel has it).
  // Isotope ratios: HITRAN92 tables.3
  species_data.push_back
    ( SpeciesRecord
      ( NAME("CO2"),
	DEGFR(3),
	ISOTOPES
	(//   Name,	Abundance,	Mass,	MY-tag, HY-tag, JPL-tag
	 //		|		|	|	|	|
	 REC( "626"	,0.98420	,44.	,21	,21	,TAGS() ),
	 REC( "636"	,0.01106	,45.	,22	,22	,TAGS() ),
	 REC( "628"	,0.0039471	,46.	,23	,23	,TAGS() ),
	 REC( "627"	,0.000734	,45.	,24	,24	,TAGS() ),
	 REC( "638"	,0.00004434	,47.	,25	,25	,TAGS() ),
	 REC( "637"	,0.00000825	,46.	,26	,26	,TAGS() ),
	 REC( "828"	,0.0000039573	,48.	,27	,27	,TAGS() ),
	 REC( "728"	,0.00000147	,47.	,28	,28	,TAGS() )
	 ) ) );
  

  // O3
  // Source: forward_4_96, glob_def.c
  // FIXME: I just copied the number without checking them. This should be checked carefully.
  //	    Do all these HITRAN tags really exist?
  //	    Are JPL tags missing?
  species_data.push_back
    ( SpeciesRecord
      ( NAME("O3"),
	DEGFR(3),
	ISOTOPES
	(//  Name,	Abundance,	Mass,	MY-tag, HY-tag, JPL-tag
	 //		|		|	|	|	|
	 REC("666"	,0.9929		,48.	,31	,31	,TAGS(48004)),
	 REC("668"	,4.0E-3		,50.	,32	,32	,TAGS(50004)),
	 REC("686"	,2.0E-3		,50.	,33	,33	,TAGS(50003)),
	 REC("667"	,7.4E-4		,49.	,34	,34	,TAGS(49002)),
	 REC("676"	,3.7E-4		,49.	,35	,35	,TAGS(49001))
	 ) ) );


}


void define_species_map()
{
}
