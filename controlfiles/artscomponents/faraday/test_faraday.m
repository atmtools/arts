% A small file to check the total Faraday rotation calculated in
% TestFaradayRotation.arts 

% 2013-02-07, Patrick Eriksson

function test_faraday

% Uncomment some lines to store the data that are required  

% Read data stored by ARTS test file  
f = xmlLoad( 'f.xml' );
n = xmlLoad( 'ne_field.xml' );
b = xmlLoad( 'bw_field.xml' );
z = xmlLoad( 'z_field.xml' ); 
r = xmlLoad( 'farrot_totalREFERENCE.xml' );

% Use integral expression shown in AUG + convert to degrees and change sign
% as trapz goes from 0 to TOA while photons travel in the opposite direction
% 
r0 = -( 180/pi * 23648 * trapz(z,b.*n) ) ./ (f.*f);

plot(f/1e9,r(1:4:end)-r0)
xlabel( 'Frequency [GHz]' )
ylabel( 'Difference [deg]' )