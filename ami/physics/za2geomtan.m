%------------------------------------------------------------------------
% NAME:    za2geomtan
%
%          Converts a vector of zenith angles to geometrical tangent altitudes.
% RETURN:  z_tan     tangent altitudes corresponding to za.
% IN:      r_geoid   Geoid radius.
%          z_plat    Platform altitude (above geoid).
%          za        zenith angles. 
%------------------------------------------------------------------------
function z_tan =  za2geomtan(r_geoid,z_plat,za);
%------------------------------------------------------------------------

DEG2RAD        = 0.01745329251994;
z_tan = (r_geoid+z_plat)*sin(DEG2RAD*za) - r_geoid;


