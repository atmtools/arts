%------------------------------------------------------------------------
% NAME:     read_mytran
%
%           Writes a Mytran line file into a structure array. 
%
%
% FORMAT:   write_mytran( filename, L )
%
% OUT:       a file
% IN:       filename   Name on file to write the data.
%           L          structure in MYTRAN format (see read_mytran,
%           for more information)
%------------------------------------------------------------------------

% HISTORY: 2003.03.28 Created by Carmen Verdes

function write_mytran( filename, L )

%= Open the file for writing
%
fid=fopen(filename,'w');

for i=1:length(L)
      
      s=[L{i}.name, L{i}.f, L{i}.df, L{i}.i0, ...
	 L{i}.agam, L{i}.sgam, L{i}.elow, L{i}.nair, ...
	 L{i}.nself, L{i}.t_gam, L{i}.psf, L{i}.V1, ...
	 L{i}.V2, L{i}.Q1,  L{i}.Q2, L{i}.di0, L{i}.dagam ...
	 , L{i}.dnair]; 
      fprintf(fid, '%s\n',  s);
end
%

%=== Close the file
fclose( fid );





function source = what_source( a )

if a == 'J'
  source = 'JPL';

elseif a == 'H'
  source = 'HITRAN';

elseif a == 'B'
  source = 'Best guess from HITRAN';

elseif a == 'T'
  source = 'Isotope default value';

elseif a == 'L'
  source = 'From default files';

else

  error(sprintf('Unknown source code (%s).',a));

end
  



function e = hitran_ecode( d )

  switch d

    case {0,1,2,3,9}
      e = -1;

    case 4
      e = 20;

    case 5
      e = 10;

    case 6
      e = 5;

    case 7
      e = 2;

    case 8
      e = 1;

    otherwise
      error(sprintf('Unknown HITRAN error code (%d)',d));
  end

return
