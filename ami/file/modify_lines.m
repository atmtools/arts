%------------------------------------------------------------------------
% NAME:     modify_lines
%
%           Modifies the intensity and broadening parameters in a line
%	    structure array.
%
%           This function gives the possibility to modify some of the
%           spectroscopic parameters. The modification is given as a scaling
%           factor. For example, to increase the value with 10%, the factor
%	    is 1.1. The input arguments starting with I0 are such scaling
%	    factors. Setting these scaling factors to be empty is equal to
%           no scaling at all (i.e. a factor of 1).
%           
%           The transitions to modify are selected by the arguments MOLEC
%	    and FREQ. Setting these arguments to be empty works as wildcards.
%           That is, setting FREQ=[] will result in that all transition
%           frequencies will be considered. If a frequency is given, it must
%           match exactly the centre frequency of the transition.
%
%           The name of the molecule will only be compared to MOLEC for the
%           length of MOLEC. This means that if only the molecule name is 
%           given (and the isotope numbers are left out), all isotopomers
%           will be considered. For example, 
%              L = modify_lines( L, 'H2O', [], 1.02 );
%           will increase the line intensity with 2% for all water 
%           transitions. 
%
% FORMAT:   [ L, nmodified ] = modify_lines( L, molec, freq, i0, agam, 
%                                                      sgam, nair, nself )
%
% OUT:      L          Structure array with line data.
%           nmodified  Number of modified transitions (lines).
% IN:       L          Structure array with line data.
%           molec      String with molecule name.
%	    freq       Frequency of transition to modify.
%           i0         Factor for line intensity ([] is allowed).
% OPTIONAL: agam       Factor for air broadened width ([] is allowed).
%           sgam       Factor for self broadened width ([] is allowed).
%           nair       Factor for temperature exponent of air broadened 
%	               width ([] is allowed).
%           nself      Factor for temperature exponent of self broadened 
%	               width ([] is allowed).
%------------------------------------------------------------------------

% HISTORY: 2002-09-06  Created by Patrick Eriksson


function [L,nmodified] = modify_lines(L,molec,freq,i0,agam,sgam,nair,nself)


if nargin < 8,   nself = [];   end
if nargin < 7,   nair  = [];   end
if nargin < 6,   sgam  = [];   end
if nargin < 5,   agam  = [];   end
if nargin < 4
  error('The function requires at least 4 input arguments.');
end
if ~( isempty( molec )  |  ischar( molec ) )
  error('The argument MOLEC must be empty or be a string.');
end


lmolec = length( molec );

do_i0    = ~isempty( i0 );
do_agam  = ~isempty( agam );
do_sgam  = ~isempty( sgam );
do_nair  = ~isempty( nair );
do_nself = ~isempty( nself );


if ~( do_i0 | do_agam | do_sgam | do_nair | do_nself )
  error('No parameter value is set to be modified.');
end


nmodified = 0;


for il = 1 : length(L)

  if isempty( molec )  |  strncmp( L{il}.name, molec, lmolec )

    if isempty( freq )  |  L{il}.f == freq

      nmodified = nmodified + 1;

      if do_i0
        L{il}.i0 = L{il}.i0 * i0;
      end

      if do_agam
        L{il}.agam = L{il}.agam * agam;
      end

      if do_sgam
        L{il}.sgam = L{il}.sgam * sgam;
      end

      if do_nair
        L{il}.nair = L{il}.nair * nair;
      end

      if do_nself
        L{il}.nself = L{il}.nself * nself;
      end

    end
  end
end