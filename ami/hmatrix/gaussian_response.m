%------------------------------------------------------------------------
% NAME:    gaussian_response
%
%          Returns a gaussian response curve.
%
%          The total width and the inverse spacing is given in units of
%          the "standard deviation" corresponding to the given FWHM.
%
%          The response curve is not normalised.
%
% FORMAT:  [x,r] = gaussian_response( fwhm, width_si, n_per_si )
%
% RETURN:  x          Abscissa (units as for WIDTH).
%          r          Response at X.
% IN:      fwhm       Full width at half maximum (in arbitrary unit).
%          width_si   The total width of the response curve in units of
%                     the "standard deviation".
%          n_per_si   Number of definition points of the response per 
%                     "standard deviation".
%------------------------------------------------------------------------

% HISTORY: 2002.01.04  Created by Patrick Eriksson. 


function [x,r] = gaussian_response( fwhm, width_si, n_per_si )


%=== Calculate the "standard deviation" corresponding to FWHM 
%
si     = fwhm/sqrt(2*log(2))/2;


%=== Create abscissa
%
width = width_si*si/2;
%
x = linspace( -width, width, round(width_si*n_per_si)+1 )';


%=== Calculate response
%
r = exp(-1*x.^2/(2*si*si))/si/sqrt(2*pi);


