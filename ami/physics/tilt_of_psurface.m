%-----------------------------------------------------------------------------
% NAME:     tilt_of_psurface
%
%           Calculates the angular tilt of a pressure surface.
%
%           The tilt is calculated numerically, and not by an analytical
%           expression.
%
%           Note that the tilt value is a local value. The tilt for a constant
%           C value, is different for different radiuses.
%
% FORMAT:   tilt = tilt_of_psurface( r, c )
%
% OUT:      tilt   The tilt of the surface in degrees.
% IN:       r      The local radius of the surface.
%           c      The slope of the surface in m/degree. 
%                  Increasing radius with incresing latitides is defined as
%                  a positive slope.
%-----------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


function tilt = tilt_of_psurface( r, c )


dalpha = 1e-2; 


tilt   = za_geom2other_point( r, r, dalpha ) - ...
                                 za_geom2other_point( r, r+c*dalpha, dalpha );

