%=============================================================
% n = vmr2nd(vmr,p,t)
%
% Calculates number density values from VMR-values.
%
% Input:	vmr	vector(s) of VMR-values [-]
%		p	vector of pressure values [Pa]
%		t	vector of temperature values [K]
%
% Output:	n	vector of number density values [1/m3]
%
% Patrick Eriksson 1993
% PE 2001-12-12, Adapted to AMI
%=============================================================

function n = vmr2nd(vmr,p,t)

global BOLTZMAN_CONST

vmr	= vec2col(vmr);
p	= vec2col(p);
t	= vec2col(t);

[rows,cols] = size(vmr);

if cols == 1
   n	= (vmr.*p./(BOLTZMAN_CONST*t)); 
else
   n	= (vmr.*(p*ones(1,cols))./(BOLTZMAN_CONST*(t*ones(1,cols)))); 
end
