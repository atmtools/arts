%------------------------------------------------------------------------
% NAME:    h_fix_ys
%
%          Returns the frequency and zenith angles vectors for the cases
%          where the frequency grid is the same for all angles.
%          This shall always be the case for the absorption and sensor
%          parts, but not always true after data reduction.
%
% FORMAT:  [f_y,za_y] = h_fix_ys(f_sensor,za_sensor)
%
% RETURN:  f_y          frequency vector        
%          za_y         zenith angle vector
% IN:      f_sensor     sensor frequencies
%          za_sensor    sensor zenith angles
%------------------------------------------------------------------------

% HISTORY: 00.08.25  Created by Patrick Eriksson. 


function [f_y,za_y] = h_fix_ys(f_sensor,za_sensor)


%=== Vector lenghts
nf        = length(f_sensor);
nza       = length(za_sensor);


%=== Allocate F_Y and ZA_Y
f_y       = zeros(nf*nza,1);
za_y      = zeros(nf*nza,1);


%=== Fill F_Y and ZA_Y 
% (This is not the most efficient Matlab solution, but it fits a future
%  implementation in ARTS.)
for i = 1:nza
  for j = 1:nf
    ind       = (i-1)*nf+j;
    f_y(ind)  = f_sensor(j);
    za_y(ind) = za_sensor(i); 
  end
end

