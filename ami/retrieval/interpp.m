%------------------------------------------------------------------------
% NAME:    interpp
%
%          Makes log-linear interpolation of pressure profiles.
%          The profiles are assumed to be constant outside the end points.
%
% FORMAT:  X = interpp(pp,Xp,p)
%
% RETURN:  X          interpolated profiles
% IN:      pp         original pressure levels
%          Xp         original profiles
%          p          new pressure levels 
%------------------------------------------------------------------------

% HISTORY: 1999.11.02  Created by Patrick Eriksson.
%          2000.01.04  Moved from Norns to AMI


function X = interpp(pp,Xp,p)

pp	= vec2col(pp);
np	= length(pp);

X	= interp1([30;log(pp);-10],[Xp(1,:);Xp;Xp(np,:)],log(p));

