%------------------------------------------------------------------------
% NAME:     read_mytran
%
%           Reads a Mytran line file into a structure array. 
%
%
% FORMAT:  [LM, LA]   = read_mytran( filename [, flims ] )
%
% OUT:      LM         Structure array for MYTRAN format. This is
%                      necessary for the case when the mitran format is
%                      request (e.g. Spectro study)
%           LA         Structure array for ARTS format
%
% IN:       filename   Name on file with Mytran line data.
% OPTIONAL: flims      Frequency limits, as [f_low,f_high], for transitions
%                      to consider.
%------------------------------------------------------------------------

% HISTORY: 2003.03.28  Created by Carmen Verdes

function [LM, LA] = read_mytran( filename, flims )


SPEED_OF_LIGHT = 2.99792458e8; 
TORR2PA        = 133.3227;
PLANCK_CONST   = 6.626180e-34;

%=== Default values
%
if nargin < 2
  flims = [0 Inf];
end


%=== Open the file     
fid = fopen( filename, 'r' );
if fid < 0
  error(sprintf('The file %s could not be opened.',filename));
end


%=== Read and store
%
i = 0;
LM = {};
LA = {};
%
s = fgets( fid );
%
while ischar( s )
  %
  if s(1) ~= '#' & length(s) > 1
    %
   

    f = sscanf( s(4:16), '%f' ) * 1e6; 

 
    if f > flims(2)
      return
    end


    if f >= flims(1)

      i = i + 1;

      %=== Reading 
  
      %= Molecule + isotope 
      %
      LM{i}.name = sscanf( s(1:3), '%c'  );
      %
      LA{i}.name = my2arts(( LM{i}.name));
      
      %= Centre frequency [Hz]
      %
      LM{i}.f = sscanf( s(4:16), '%c'  );
      %
      LA{i}.f = str2num( LM{i}.f ) * 1e6;
      
      % Accuracy of frequency
      %
      LM{i}.df = sscanf( s(17:24), '%c'  );
      %
      LA{i}.df = str2num( LM{i}.df ) * 1e6;
      
      % Intensity in cm-1/(mole cm-2) at 296K
      %
      LM{i}.i0 = sscanf( s(25:34), '%c'  );
      %      
      hi2arts = 1e-2 * SPEED_OF_LIGHT;
      %
      LA{i}.i0 = str2num( LM{i}.i0) * hi2arts;
      
      % Reference temperature for i0
      % 
      LA{i}.t_i0 = 296; 
      
      % Agam in MHz/Torr at T0
      %
      LM{i}.agam = sscanf( s(35:39), '%c'  );
      %
      LA{i}.agam =str2num( LM{i}.agam )* 1e6/ TORR2PA;
      
      % Sgam in MHz/Torr at T0
      %
      LM{i}.sgam = sscanf( s(40:44), '%c'  );
       %
      LA{i}.sgam =str2num( LM{i}.sgam )* 1e6/ TORR2PA;
      
      % Lower state energy in wavenumbers (cm-1)
      %
      LM{i}.elow = sscanf( s(45:54), '%c'  );  
      %
      wavenumber_to_joule = PLANCK_CONST * SPEED_OF_LIGHT * 1e2;
      LA{i}.elow = wavenumber_to_joule * str2num( LM{i}.elow );
      
      % Nair
      %
      LM{i}.nair = sscanf( s(55:58), '%c'  );
      %
      LA{i}.nair = str2num( LM{i}.nair ); 
      
      % Nself
      %
      LM{i}.nself = sscanf( s(59:62), '%c'  );
      %
      LA{i}.nself = str2num( LM{i}.nself ); 
      
      
      % T0, reference temperature for Agam and Sgam
      %
      LM{i}.t_gam = sscanf( s(63:69), '%c'  ); 
      % 
      LA{i}.t_gam = str2num( LM{i}.t_gam );
      
      %= Pressure shift parameter [MHz/Torr]
      %
      LM{i}.psf = sscanf( s(70:78), '%c' );
      %
      LA{i}.psf =  str2num( LM{i}.psf ) * 1e6 /TORR2PA;
      
      % Upper state global quanta index
      %
      LM{i}.V1 = sscanf( s(79:81), '%c' );
      %
      LA{i}.qcode = str2num( LM{i}.V1 );
      
      % Lower state global quanta index
      %
      LM{i}.V2 = sscanf( s(82:84), '%c' );
      %
      % this field does not exist in ARTS
      
      % Upper state local quanta
      %
      LM{i}.Q1 = sscanf( s(85:93), '%c' );
      %
      LA{i}.qlower = LM{i}.Q1;
  
      % Upper state local quanta
      %
      LM{i}.Q2 = sscanf( s(94:102), '%c' );
      %
      LA{i}.qupper = LM{i}.Q2;
      
      % Accuracy index for I0
      %
      LM{i}.di0 = sscanf( s(103), '%c' );
      %
      LA{i}.di0 = convMytranIER (str2num( LM{i}.di0 ) );
      
      % Accuracy index for agam
      %
      LM{i}.dagam = sscanf( s(104), '%c' );
      %
      LA{i}.dagam = convMytranIER (str2num( LM{i}.dagam ) );
      
      % Accuracy index for nair
      %
      LM{i}.dnair= sscanf( s(105), '%c' );
      %
      LA{i}.dnair = convMytranIER (str2num( LM{i}.dnair ) );
      
      % Accuracy index for sgam
      %
      % This does not exist in MYTRAN
      %
      LA{i}.dsgam = -1;
      
      % Accuracy index for nself
      %
      % This does not exist in MYTRAN
      %
      LA{i}.dnself = -1;
            
      %=== Number of aux parameters [-]
      %
      LA{i}.n_aux = 0;
      
      % accuracy for pressure shift
      %
      LA{i}.dpsf = -1;
           
    end
  end
  %
  s = fgets( fid );
  %
end


%=== Close the file
fclose( fid );


%----------------------------------------------------------------
function  da = convMytranIER (dm)
%convert MYTRAN index for intensity and halfwidth accuracy to ARTS
%units (relative units).

switch ( dm )
 case 0 
  da = 20;
 case 1
  da = 100; 
 case 2
  da = 50;
 case 3
  da = 30; 
 case 4
  da = 20;
 case 5
  da = 10;
 case 6
  da =5;
 case 7
  da = 2;
 case 8
  da = 1;
 case 9
  da = 0.5;
 otherwise
  display ('nonexistent error code');
end

da = da / 100;

