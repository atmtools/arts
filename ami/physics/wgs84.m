%------------------------------------------------------------------------
% NAME:    wgs84
%
%          Calculates a radius of the WGS-84 reference geiod ellipsoid
%          suitable to use as the geoid radius in a atmospheric forward
%          model.
%
%          For atmospheres with a spherical symmetry (DIM=1), the 
%          curvature radius at the given latitude(s) is returned. This
%          radius is the circle radius matching the local shape of the
%          reference ellipsoid. The equations are taken from the book
%          by Rodgers (Eqs. 9.27-9.29).
%
%	   For atmospherers not assuming a spherical Earth (DIM=2 and 
%          DIM=3), the ellipsoid radius is returned.
%
%          The OBSDIR argument is only a valid input argument for DIM=1.
%
%          The latitude argument can be a vector, while the observation
%	   direction can only be a scalar.
%
% FORMAT:  r = wgs84( dim, lat, obsdir )
%
% RETURN:  r        Geoid radius.
% IN:      dim      Dimension of simulated atmopshere, 1-3.
%          lat      Geocentric latitude
%          obsdir   Azimuth angle to the meridian plane. 
%------------------------------------------------------------------------

% HISTORY: 2002-02-08  Created by Patrick Eriksson. 
%          2002-03-06  Added the DIM argument (PE).


function r = wgs84( dim, lat, obsdir )


global DEG2RAD

rq = 6378.138e3;
rp = 6356.752e3;

rq = rq * rq;
rp = rp * rp;

sv = sin( lat * DEG2RAD );
cv = cos( lat * DEG2RAD );

sv = sv .* sv;
cv = cv .* cv;


if dim<1 | dim>3
  error('Valid values for DIM are 1-3.');
end

if dim == 1

  if nargin < 3
    error('For DIM=1, an azimuth angle must be specified.');
  end

  v = rq*cv + rp*sv;

  rns = rq * rp * v.^-1.5;
  rew = rq * v.^-0.5;

  sv = sin( obsdir * DEG2RAD );
  cv = cos( obsdir * DEG2RAD );

  r = 1 ./ ( cv*cv./rns + sv*sv./rew );


else

  if nargin > 2
    error('To many input arguments for DIM=2 or DIM=3.');
  end

  r = sqrt( rq*rp ./ ( rq*sv + rp*cv ) );

end