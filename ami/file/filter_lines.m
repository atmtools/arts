%------------------------------------------------------------------------
% NAME:     filter_lines
%
%           Filter a line array to remove unwanted lines.
%
%           The lines are given as structure array. Field names for the 
%           structures are defined in WRITE_LINEFILE.
%
%           A single filter can be given as string, for example:
%
%              L = filter_lines( L, 'H2O-PWR98' );
%
%           If several filters shall be applied, the filters are specified
%           as a cell array with strings:
%
%              L = filter_lines( L, {'H2O-PWR98','O2-PWR98'} );
%
%           If the filters take any arguments, these are given as optional
%           cell arrays. A maximum of two filter arguments are accpted:
%
%              L = filter_lines( L, 'i0_threshold', {'O3-666'}, {1e-17} );
%
%           Dummy values are accepted for filters not taken any arguments.
%           Note that the filter argum,ents always are cell arrays, even if
%           only a single filter is applied.
%
%           The filters are implemented inside the function. See the file
%           for instructions on how to add more filters.
%
%           The following filters are implemented:
%
%           i0_threshold: 
%           Removes transitions of an isotope with an intensity below the 
%           given threshold. The isotope name shall be given as first extra 
%           argument, and the intensity threshold as second. For  example: 
%              filter_lines( L, 'i0_threshold', {'O3-666'}, {1e-17} );
%
%           H2O-PWR98: 
%           Removes all H2O-161 transitions that are included in the water 
%           vapour absorption model by Rosenkrantz 1998.
%
% FORMAT:   Lout = filter_lines( Lin, filter )
%
% OUT:      Lout     Output structure array with line data.
% IN:       Lin      Input structure array with line data.
%           filter   String or cell array with filter names.
%------------------------------------------------------------------------

% HISTORY: 2002.01.03  Created by Patrick Eriksson

% Filters are implemented below the dashed the line, where each filter
% is a seperate case statement.
%
% The filter can be implemented in any way as long as the variable OK is set
% to 0 (zero) for lines that shall be removed. The variable OK is set to 1 as
% default.
% As described above, the filters can take two arguments. These arguments
% are passed to the filter sub-function as FILTARG1 and FILTARG2.
% Good luck. Patrick Eriksson 2002-01-03.


function Lout = filter_lines( Lin, filter, filtarg1, filtarg2 )


%=== Check type of FILTER, convert to cell if FILTER is a string
%
if ischar(filter)
  %
  filter = {filter};

elseif iscell(filter)
  %
  % Nothing to do

else
  error('FILTER must be a cell array or a string');
end


%=== Create dummy fiter arguments if not already defined
%
nfilters = length(filter);
%
if nargin < 3
  filtarg1 = cell( nfilters, 1 );
end
if nargin < 4
  filtarg2 = cell( nfilters, 1 );
end



%=== Counter for lines not removed
%
for ii = 1 : nfilters

  if ~ischar( filter{ii} )
    error(sprintf('Filter name %d is not a string',ii));
  end

  if ii > 1
    Lin = Lout;
  end

  nout = 0;
  Lout = [];
  
  for i = 1 : length(Lin)
    %
    if linefilter( Lin{i}, filter{ii}, filtarg1{ii}, filtarg2{ii} )
      %
      nout = nout + 1;
      Lout{nout} = Lin{i};
      %
    end
    %
  end
end
return


%----------------------------------------------------------------------------


function ok = linefilter( L, filter, filtarg1, filtarg2 )

%= Set OK=1 as default
%
ok = 1;

switch filter


%--- i0_treshold --------------------------------------------------------------
%
% Intensity treshold. Description of the filter in the function header.
%
%
case 'i0_threshold'
  %
  if ~ischar( filtarg1 )
    error('The first argument to i0_threshold shall be a string.');
  end
  %
  if strcmp( L.name, filtarg1 )
    %
    if ~isscalar( filtarg2 )
      error('The second argument to i0_threshold shall be a scalar.');
    end
    %
    if L.i0 < filtarg2
      ok = 0;
    end
    %
  end



%--- H2O-PWR98 --------------------------------------------------------------
%
% Water wapour absorption model by Rosenkrantz 1998.
%
% This model includes only transitions of the main isotope. 
% A match is considered to be found if the transition frequency is closer 
% than 1 MHz to the lines in the absorption model (see below). 
% The smallest frequency seperation between two H2O-161 lines below 1 THz
% in the JPL catalogue is 18 MHz.
% Note that the lines variable must be a vector (not a matrix).
%
case 'H2O-PWR98'
  %
  if strcmp( L.name, 'H2O-161' )
    %
    lines = [ 22.2350800,183.3101170,321.2256400,325.1529190,380.1973720, ...
             439.1508120,443.0182950,448.0010750,470.8889470,474.6891270, ...
             488.4911330,556.9360020,620.7008070,752.0332270,916.1715820]*1e9;
    %
    if any( abs( L.f - lines ) <= 1e6 )
      ok = 0;
    end 
    %
  end


%--- Unknown filter ---------------------------------------------------------
%
otherwise
  error(sprintf('Unknown line filter (%s)',filter));
end
return
