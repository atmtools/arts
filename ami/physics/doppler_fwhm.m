%=============================================================
% w = doppler_fwhm(f0,t,molweight)
%
% Calculates the full width at half maximum of Doppler
% broadening. Expression taken from my (PE) thesis, page 221.
%
% Input:	f0	    centre frequency of transition 
%		t	    temperature
%		molweight   molecular weight (e.g. 48 for O3)
%
% Output:	w	FWHM of Doppler broadening
%
% Patrick Eriksson 2003-09-11
%=============================================================

function w = doppler_fwhm(f0,t,molweight)

w = 7.16e-7 * f0 * sqrt( t / molweight );
