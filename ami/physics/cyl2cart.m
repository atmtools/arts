%------------------------------------------------------------------------
% NAME:     cyl2cart
%
%           Converts polar coordinates to cartesian coordinates.
%
%           The zero point for the angle in the polar coordinate system is
%           along the y-axis. 
%
% FORMAT:   [x,y] = cyl2cart( r, alpha )
%
% OUT:      x       X-coordinate in the cartesin system.
%           y       Y-coordinate in the cartesin system.
% IN:       r       The radius in the polar coordinate system.
%           alpha   Angle from the y-axis.
%------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


function [x,y] = cyl2cart( r, alpha )


global DEG2RAD

x = r * sin( DEG2RAD*alpha );
y = r * cos( DEG2RAD*alpha );

