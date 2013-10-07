% A small file to check the total Faraday rotation calculated in
% TestFaradayRotation.arts 

% 2013-02-07, Patrick Eriksson

function test_faraday

% Uncomment some lines to store the data that are required  

% Read data stored by ARTS test file  
f   = xmlLoad( 'f.xml' );
vmr = xmlLoad( 'vmr_field.xml' );
b   = xmlLoad( 'bw_field.xml' );
z   = xmlLoad( 'z_field.xml' ); 
r   = xmlLoad( 'farrot_totalREFERENCE.xml' );
%
n   = vmr(end,:)';


% Use integral expression shown in AUG + convert to degrees and change sign.
% The later as magnetic field vector and photon direction are reversed, and
% the angular part of the dot product equals -1. 
% 
r0 = -( 180/pi * 2.364797970062947e+04 * trapz(z,b.*n) ) ./ (f.*f);

figure(1)
plot(f/1e9, r(1:4:end)-r0 )
xlabel( 'Frequency [GHz]' )
ylabel( 'Difference [deg]' )

figure(2)
plot(f/1e9, r(1:4:end) )
xlabel( 'Frequency [GHz]' )
ylabel( 'Calculated rotation [deg]' )

