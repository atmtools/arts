%------------------------------------------------------------------------
% NAME:     write_linefile
%
%           Creates a line file in the ARTS format.
%
%           The line data shall be given by an array of structures. The
%           structures shall have fields corresponding to the items
%           of the line format. The field names are in lower case.
%           No unit conversions are performed, thus the data shall have 
%           the smae unit as in the line file.
%
% FORMAT:   write_linefile( filename, L [, do_quanta, do_source ])
%
% OUT:      -
% IN:       filename   Name on file to create.
%           L          Structure array with line data (see above)
% OPTIONAL: do_quanta  Flag for writing of quantum data. Default is 1. 
%           do_source  Flag for writing of source information. Default is 1.
%------------------------------------------------------------------------

% HISTORY: 2001.03.27  Created by Patrick Eriksson


function write_linefile( filename, L, do_quanta, do_source )

if ~exist('do_quanta'),   do_quanta = 1;   end
if ~exist('do_source'),   do_source = 1;   end


%=== Open the file     
fid = fopen( filename, 'w' );
if fid < 0
  error(sprintf('The file %s could not be opened.',filename));
end


%=== Write data
%
for i = 1:length(L)

  fprintf(fid, '@ %s',   L{i}.name );
  fprintf(fid, ' %.9e',  L{i}.f );
  if L{i}.psf == 0
    fprintf(fid, ' 0');
  else  
    fprintf(fid, ' %.6e',  L{i}.psf );
  end
  fprintf(fid, ' %.6e',  L{i}.i0 );
  fprintf(fid, ' %.0f',  L{i}.t_i0 );
  fprintf(fid, ' %.4f',  L{i}.elow );
%  fprintf(fid, '\n   ');
  fprintf(fid, ' %.2f',  L{i}.agam );
  fprintf(fid, ' %.2f',  L{i}.sgam );
  fprintf(fid, ' %.2f',  L{i}.nair );
  fprintf(fid, ' %.2f',  L{i}.nself );
  fprintf(fid, ' %.0f',  L{i}.t_gam );
  fprintf(fid, ' %d',    L{i}.n_aux );
  for j = 1:L{i}.n_aux
    fprintf(fid, ' %.6e',  eval(['L{i}.aux',int2str(j)]) );
  end
%  fprintf(fid, '\n   ');
  fprintf(fid, ' %.3e',  L{i}.df );
  fprintf(fid, ' %.1f',  L{i}.di0 );
  fprintf(fid, ' %.1f',  L{i}.dagam );
  fprintf(fid, ' %.1f',  L{i}.dsgam );
  fprintf(fid, ' %.1f',  L{i}.dnair );
  fprintf(fid, ' %.1f',  L{i}.dnself );
  fprintf(fid, ' %.1f',  L{i}.dpsf );

  if do_quanta
%    fprintf(fid, '\n   ');
    fprintf(fid, ' "%s"',  L{i}.qcode );
    fprintf(fid, ' "%s"',  L{i}.qlower );
    fprintf(fid, ' "%s"',  L{i}.qupper );

    if do_source
%      fprintf(fid, '\n   ');
      fprintf(fid, ' "%s"',  L{i}.if );
      fprintf(fid, ' "%s"',  L{i}.ii0 );
      fprintf(fid, ' "%s"',  L{i}.ilw );
      fprintf(fid, ' "%s"',  L{i}.ipsf );
      fprintf(fid, ' "%s"',  L{i}.iaux );
    end
  end

  fprintf(fid, '\n');


end


%=== Close the file
fclose( fid );
