%-----------------------------------------------------------------------------
% NAME:     zpot2zgeom
%
%	    Converts potential altitudes to geomtrical altitudes.
%
% FORMAT:   z = zpot2zgeom( z, lat )
%
% OUT:      z     Geometrical altitude.
% IN:       z     Potential altitude.
%           lat   Latitude.
%-----------------------------------------------------------------------------

% HISTORY: 2002-08-14  Created by Patrick Eriksson


function z = zpot2zgeom(z,lat)


r = wgs84( 3, lat );


z = r * z ./ ( r - z ); 