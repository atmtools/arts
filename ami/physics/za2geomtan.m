%------------------------------------------------------------------------
% NAME:    za2geomtan
%
%          Converts a vector of zenith angles to geometrical tangent
%          altitudes 
%
% FORMAT:  z_tan = geomtan2za(r_geoid,z_plat,za)
%
% RETURN:  z_tan     Geometrical tangent altitudes. 
% IN:      r_geoid   Geoid radius.
%          z_plat    Platform altitude (above geoid).
%          za        Zenith angles corresponding to z_tan.
%------------------------------------------------------------------------

% HISTORY: 2001.10.15  Created by Carlos Jimenez 


function z_tan = za2geomtan(r_geoid,z_plat,za)


z_tan = (r_geoid+z_plat).*sin(pi*(180-za)/180)-r_geoid;

