%=============================================================
% vmr = nd2vmr(n,p,t)
%
% Calculates VMR-values from number density values.
%
% Input:	n	vector of number density values [1/m3]
%		p	vector of pressure values [Pa]
%		t	vector of temperature values [K]
%
% Output:	vmr	vector of VMR-values [ppm]
%
% Patrick Eriksson 1993
% PE 2001-12-12, Adapted to AMI
%=============================================================

function vmr = nd2vmr(n,p,t)

global BOLTZMAN_CONST

vmr	= n*BOLTZMAN_CONST.*t./p;
