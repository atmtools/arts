%-----------------------------------------------------------------------------
% NAME:     cira86
%
%           Extracts data from the COSPAR International Reference Atmosphere,
%           CIRA-86.
%
%	    The function requieres that arts has been configured with 
%           arts-data. If not, there should have been a warning when the AMI
%	    init script was run.
%
%           The data are read from the files sht.dat, nht.dat, shz.dat and
%           nhz.dat and are interpolated to the given points. The interpolation
%           is performed with interp3, where the linear option is used. That
%           is, 3D linear interpolation. Note that fractional month numbers
%           are allowed. Allowed month values range from 0 to 13.
%
%	    The interpolated data are returned as tensor of order 3, where 
%	    the index are (using temperature as example):
%
%              T(month,latitude,pressure)
%
%           Data for latitudes <= -70 and low altitudes are missing (there is 
%           only ice there). The temperature is set to 250 K for the missing
%           data points, and some approximative altitudes are calculated, to
%           make the data set complete.
%
%           Note that the function returns geometrical altitudes, not potential
%           altitudes as original CIRA data. Note further that the CIRA 
%	    altitudes do not fulfill hydrostatic eq. perfectly.
%
%           The input vectors must be sorted.
%
% FORMAT:   [T,Z] = cira86(months,lats,p)
%
% OUT:      T       Temperature          [month,lat,p]
%           Z       Geometrical altitude [month,lat,p]
% IN:       months  A vector with month numbers (1=Jan etc.).
%           lats    A vector with latitudes.
%           p       A vector with pressure values [Pa]
%-----------------------------------------------------------------------------

% HISTORY: 2002-08-14  Created by Patrick Eriksson


function [T,Z] = cira86(months,lats,p)


%=== Check input 
if any( months < 0 )  |  any( months > 13 )
  error('Allowed range for month numbers is 1 - 12.');
end
if any( lats < -90 )  |  any( lats > 90 )
  error('Allowed range for latitudes is -90 - 90.');
end
if any( lats < -90 )  |  any( lats > 90 )
  error('Allowed range for latitudes is -90 - 90.');
end
if any( p < 2.54e-3 )  |  any( p > 1.013e5 )
  error('Allowed range for pressures is 2.54e-3 - 101300 Pa.');
end


%=== Determine the path to the CIRA folder in arts-data.
%
if ~exist('arts_data_path','file');
  error( ['It seems that init.m could not find arts-data, that is ',...
                                               'required for this function'] );
end
%
cira_path = fullfile( arts_data_path, 'atmosphere', 'cira86' );




%=== Latitude and month grids for read CIRA data (pressures are read from
%=== one of the files)
%
clats = [-90 -80:5:80 90]';
cmonths = [0:13]';



%=== Read southern hemisphere temperatures
%
fid = cira_open( cira_path, 'sht.dat' );
[cp,SH] = cira_read( fid, 't' );
fclose(fid);


%=== Read northern hemisphere temperatures
%
fid = cira_open( cira_path, 'nht.dat' );
[cp,NH] = cira_read( fid, 't' );
fclose(fid);


%=== Combine the temperatures for SH and NH
%
TT = cira_combine( SH, NH );


%=== Fix bad values at the lowest two altitudes and latitudes <= -70
%
TT = cira_fix_bad_t( TT );


%=== Interpolate
%
T = interp3(clats,cmonths,cp,TT,lats,months,p);



%=== Do the same for geopotential altitudes
%
if nargout > 1

  fid = cira_open( cira_path, 'shz.dat' );
  [cp,SH] = cira_read( fid, 'z' );
  fclose(fid);

  fid = cira_open( cira_path, 'nhz.dat' );
  [cp,NH] = cira_read( fid, 'z' );
  fclose(fid);

  TT = cira_combine( SH, NH );

  TT = cira_fix_bad_z( TT );

  Z = interp3(clats,cmonths,cp,TT,lats,months,p);

  % Convert to geometrical altitudes
  for i = 1:length(lats);
    Z(:,i,:) = zpot2zgeom( Z(:,i,:), lats(i) );
  end


end


%----------------------------------------------------------------------------

function fid = cira_open( folderpath, name )
  fid = fopen( fullfile(folderpath,name), 'r' );
  if fid < 0
    error(sprintf('The file %s could not be opened.',name));
  end
return



function [cp,A] = cira_read( fid, t_or_z )

  if t_or_z == 't'
    n = 20;
  elseif t_or_z == 'z'
    n = 19;
  else
    error('Unknown letter code');
  end

  % Allocate a tensor for data to read
  B = zeros(12,n,71);

  % Loop months and read
  for i = 1:12

    % Read 6 lines with text
    for j = 1:6
      fgetl(fid);
    end 

    % Read the data for the month
    B(i,:,:) = fscanf( fid, '%f', [n 71] ); 

    % Read newline character
    fgetl(fid);

  end

  cp = squeeze( B(1,2,:) ) * 100;
  A  = B( :, n + (-16:0), : );
  
return



function A = cira_combine( SH, NH )

  A = zeros(14,35,71);

  % Data for +-80 are copied to be valid also for +-90
  % Data for month 1 is copied to month 13, and data for 12 to 0.

  A(2:13,1,:)     = SH(:,1,:);
  A(2:13,35,:)    = NH(:,17,:);
  A(2:13,2:18,:)  = SH;
  A(2:13,19:34,:) = NH(:,2:17,:);

  A(1,:,:)        = A(13,:,:);
  A(14,:,:)       = A(2,:,:);
return



function T = cira_fix_bad_t( T )

  % Set the bad values (=999.9) to 250 K
  T( find( T == 999.9 ) ) = 250;

return



function Z = cira_fix_bad_z( Z )

  % Assume 1900 m between the pressure levels and go downwards
  for m = 1:12
    for l = 1:4
      for p = 70:71
        Z(m,l,p) = Z(m,l,p-1) - 1900;
      end
    end
  end

return