%------------------------------------------------------------------------
% NAME:     combine_lines
%
%           Combines "identical" transitions to a single transition.
%
%           To consider transitions to be "identical", their frequencies
%           must be inside DF, and the the lower energy and the air 
%           broadening line width must coincide inside the specified 
%           relative limits.
%
%           The intensity of the new transition is the sum of the 
%           intensity of the combined transitions.
%
%           If exist, the fields QLOWER, QUPPER and II0 are changed to
%           flag that the new transition is a combination of lines.
%
% FORMAT:   Lout = combine_lines( Lin, df [, delow, dagam ] )
%
% OUT:      Lout    Structure array with line data (see WRITE_LINEFILE).
% IN:       Lin     Structure array with line data (see WRITE_LINEFILE).
%           df      Maximum difference in transition frequency for 
%                   considering to combine the transitions.
% OPTIONAL: delow   Maximum allowed relative difference in lower energy
%                   to combine transitions. Default is 1e-6;
% OPTIONAL: dagam   Maximum allowed relative difference in air broadening
%                   line width to combine transitions. Default is 1e-6;
%------------------------------------------------------------------------

% HISTORY: 2002.01.15  Created by Patrick Eriksson


function Lout = combine_lines( Lin, df, delow, dagam )


%=== Handle optional arguments
%
if nargin < 3
  delow = 1e-6;
end
if nargin < 4
  dagam = 1e-6;
end


%=== Allocate Lout
%
Lout = Lin;
%
nout = 1;


%=== Loop lines
%
for i = 2 : length(Lin)

  if abs(Lin{i}.f-Lout{nout}.f) < df & ...
     abs((Lin{i}.elow-Lout{nout}.elow)/Lout{nout}.elow) < delow & ...
     abs((Lin{i}.agam-Lout{nout}.agam)/Lout{nout}.agam) < dagam 

    Lout{nout}.i0 = Lout{nout}.i0 + Lin{i}.i0;

    if isfield( Lin{i}, 'qlower' )
      Lout{nout}.qlower = 'combined';
      Lout{nout}.qupper = 'combined';
    end
    if isfield( Lin{i}, 'ii0' )
      Lout{nout}.ii0    = sprintf('%s summed',Lin{i}.ii0);
    end

  else

    nout = nout + 1;
    %
    Lout{nout} = Lin{i};

  end

end


if nout < length(Lin)

  Lout = Lout(1:nout);

end