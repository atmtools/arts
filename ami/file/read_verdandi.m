%------------------------------------------------------------------------
% NAME:     read_verdandi
%
%           Reads a Verdandi line file into a structure array. 
%
%           The line file can have comments lines starting with # and
%           empty lines just including a line break.
%
%           The Verdandi format is, for example, found in my Ph.D. thesis
%           on page 229.
%
%           Each structure in the array have a field for each item in the
%           ARTS format. The units are also exactly as for the ARTS format.
%           See further write_linefile in AMI.
%
%           The HITRAN codes are used for the transition intensity and
%           line widths. Codes 0-3 (where 3 means >20%) and 0 are translated
%           to -1 to indicate that no error estimates exist. 
%
%           The optional argument FLIMS can be used to only pick out 
%           transitions inside a frequency range. This option is best used
%           when reading a file with frequency sorted transitions as the
%	    reading is stoped as soon as a transition is above the upper
%           frequency limit.
%
% FORMAT:   L = read_verdandi( filename [, flims ] )
%
% OUT:      L          Structure array (see above)
% IN:       filename   Name on file with Verdandi line data.
% OPTIONAL: flims      Frequency limits, as [f_low,f_high], for transitions
%                      to consider.
%------------------------------------------------------------------------

% HISTORY: 2001.03.27  Created by Patrick Eriksson

function L = read_verdandi( filename, flims )


%=== Default values
%
if nargin < 2
  flims = [0 Inf];
end


%=== Conversion from JPL/Verdandi intensity unit to m2/Hz
%
ijpl2m2hz = 1e-12;


%=== Conversion factor from MHz/Torr to Hz/Pa
%
mhztorr2hzpa = 1e6* 760 / 1.013e5;


%=== Conversion factor from cm-1 to J
%
cminv2j = 29.97925e9 * 6.6262e-34;



%=== Open the file     
fid = fopen( filename, 'r' );
if fid < 0
  error(sprintf('The file %s could not be opened.',filename));
end


%=== Read and store
%
i = 0;
L = {};
%
s = fgets( fid );
%
while ischar( s )
  %
  if s(1) ~= '#' & length(s) > 1
    %
    if length(s) < 110
      error(sprintf('Line %d seems to be imcomplete.'));
    end

    f = sscanf( s(4:16), '%f' ) * 1e6; 

 
    if f > flims(2)
      return
    end


    if f >= flims(1)

      i = i + 1;

      %=== Reading and conversion
  
      %= Species name    
      %
      L{i}.name = vtag2arts( sscanf( s(1:3), '%f' ) );
  
      %= Centre frequency [Hz]
      %
      L{i}.f = f;
  
      %= Pressure shift parameter [Hz/Pa]
      %
      L{i}.psf = sscanf( s(25:32), '%f' ) * mhztorr2hzpa;  
  
      %= Line intensity [m^2/Hz]
      %
      L{i}.i0 = 10^sscanf( s(33:40), '%f' ) * ijpl2m2hz;  
  
      %= reference temperature for I0
      %
      L{i}.t_i0 = 300;
  
      %=== Lower state energy [cm-1]
      %
      L{i}.elow = sscanf( s(41:50), '%f' ) * cminv2j;  
  
      %=== Air broadened half-width [Hz/Pa]
      %
      L{i}.agam = sscanf( s(51:55), '%f' ) * mhztorr2hzpa;  
  
      %=== Self broadened half-width [Hz/Pa]
      %
      L{i}.sgam = sscanf( s(56:60), '%f' ) * mhztorr2hzpa;  
  
      %=== AGAM temperature exponent [-]
      %
      L{i}.nair = sscanf( s(61:64), '%f' );  
  
      %=== SGAM temperature exponent [-]
      %
      L{i}.nself = sscanf( s(65:68), '%f' );  
  
      %=== Reference temperature for AGAM and SGAM [K]
      %
      L{i}.t_gam = 296;  
  
      %=== Number of aux parameters [-]
      %
      L{i}.n_aux = 0;  
  
      %= Uncertianty of centre frequency [Hz]
      %
      L{i}.df = sscanf( s(17:24), '%f' ) * 1e6;  
  
      %=== Error for I0 [%]
      %
      a = s(99);  
      L{i}.di0 = -1;
  
      %=== Error for AGAM [%]
      %
      L{i}.dagam = hitran_ecode( str2num(s(100)) );
  
      %=== Error for SGAM [%]
      %
      L{i}.dsgam = hitran_ecode( str2num(s(101)) );
  
      %=== Error for NAIR [%]
      %
      L{i}.dnair = hitran_ecode( str2num(s(102)) );
  
      %=== Error for NSELF [%]
      %
      L{i}.dnself = hitran_ecode( str2num(s(103)) );
  
      %=== Error for PSF [%]
      %
      L{i}.dpsf = -1;
  
      %=== Quantum number code (string)
      %
      L{i}.qcode = s(71:74);  
  
      %=== Lower state quanta (string)
      %
      L{i}.qlower = s(75:86);  
  
      %=== Upper state quanta (string)
      %
      L{i}.qupper = s(87:98);  
  
      %=== Information source for F [string]
      %
      L{i}.if = what_source( s(104) );
  
      %=== Information source for I0 [string]
      %
      L{i}.ii0 = what_source( s(105) );
  
      %=== Information source for line width parameters [string]
      %
      L{i}.ilw = what_source( s(106) );
  
      %=== Information source for pressure shift [string]
      %
      L{i}.ipsf = what_source( 'H' );
  
      %=== Information source for aux [string]
      %
      L{i}.iaux = '-';

    end
  end
  %
  s = fgets( fid );
  %
end



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
