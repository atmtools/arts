%------------------------------------------------------------------------
% NAME:    hseCalc
% 
%          Calculates hydrostatic equilibrium resembling the ARTS
%          hse calculation
%          
% FORMAT:  z_abs = hseCalc(p_abs,t_abs,z_abs,h2o_abs,pref,zref,r_geoid,niter)
%
% RETURN:  z_abs      altitude profile fulfilling hse with t_abs and z_abs
%                     input
%
% IN:      p_abs      pressure vector
%          t_abs      temperature vector
%          z_abs      altitude vector
%          h2o_abs    H2O profile in p_abs grid
%          pref       reference pressure
%          zref       reference altitude
%          r_geoid    radius of the Earth's geoid corresponding to
%                     the geografical location of the given profiles
%          niter      number of iterations, 2 should be enough
%     
%-------------------------------------------------------------------------

% HISTORY: 01.10.07  Created by Carlos Jimenez. 


function z_abs = hseCalc(p_abs,t_abs,z_abs,h2o_abs,pref,zref,r_geoid,niter)


i        = 0;             % altitude index
g        = 0;             % gravitational acceleration
r        = 0;             % water mixing ratio in gram/gram
tv       = 0;             % virtual temperature
dz       = 0;             % step geometrical altitude
g0       = 9.81;          % gravitational acceleration

np       = length(p_abs); 
ztmp     = zeros(np,1);   % temporary storage for z_abs



for iter=1:niter

      %-- Init ztmp
      ztmp(1) = z_abs(1);
  
      %-- Calculate new altitudes (relative z_abs(1)) and store in ztmp
      for i=1:np-1
      
	%-- Calculate g 
	g  = ( g_of_z(r_geoid,g0,z_abs(i)) + ...
	       g_of_z(r_geoid,g0,z_abs(i+1)) ) / 2;
  
	%-- Calculate weight mixing ratio for water assuming constant average
	%   molecular weight of the air
	r  = (18/28.96) * (h2o_abs(i)+h2o_abs(i+1)) / 2;
  
	%--  The virtual temperature (no liquid water)
	tv = (1+0.61*r) * (t_abs(i)+t_abs(i+1)) / 2;

  
	%-- The change in vertical altitude from i to i+1 
	dz = 287.053 * (tv/g) * log( p_abs(i)/p_abs(i+1) );
	ztmp(i+1) = ztmp(i) + dz;
      
      end
  

      %-- Match the altitude of the reference point
      dz = interpp( p_abs, ztmp, pref ) - zref;

      z_abs = ztmp - dz;

end

return

%---------------------------------------------------------------------------
%
%                                 SUBFUNCTIONS
%
%---------------------------------------------------------------------------


function g = g_of_z (r_geoid,g0,z)

g = g0 * power( r_geoid/(r_geoid+z), 2 );



