%------------------------------------------------------------------------
% NAME:    hFreqSwitch
%
%          Includes in H the result of a frequency switching
%          (no folding).
%          
%          The given frequency shift is the total frequency shift
%          assuming a symmetric throw. Accordingly, an observed line with
%          center frequency F0, will appear at F0-F_SHIFT/2 and at 
%          F0+F_SHIFT/2. 
%
%          The defualt is that the line appearing at F0-F_SHIFT/2 is weighted
%          negetivaly (thus looking as an absorptionb line). This can be
%          changed be setting the optional argument RWEIGHTS.
%
%          The spacing between frequency channels must equally spaced, and
%          the frequency shift must be a multiple of this spacing.
%
%          Only observations with one viewing angle are handled.
%
% FORMAT:  [H,f_y,za_y,f_sensor] = hFreqSwitch(H,f_sensor,za_sensor,f_shift 
%                                                                 [,rweights] )
%
% RETURN:  H           H matrix after frequency shift
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          f_sensor    new frequency grid 
% IN:      H           H matrix after the mixer
%          f_sensor    input frequency grid
%          za_sensor   zenith angles
%          f_shift     full frequency shift
% OPTIONAL:rweights    reveresed the weighting between the 2 frequency throws
%------------------------------------------------------------------------

% HISTORY: 01.12.08  Created by Patrick Eriksson 


function [H,f_y,za_y,f_sensor] = hFreqSwitch(H,f_sensor,za_sensor,f_shift,rweights)


%=== Only one viewing angle is handled
%
if length(za_sensor) > 1
  error('Frequency shift is only handled for 1 viewing angle.');
end


%=== Check frequency shifts
%
if length(f_shift) ~= 1
  error('The frequency shift shall be a scalar, not a vector.');
end


%=== Frequency vector must be equally spaced
%
df = f_sensor(2) - f_sensor(1);
%
if any( diff(f_sensor) ~= df )
  error('The input frequency vector must be equally spaced.');
end


%=== Determine how many channels the shift equals
%
if abs(rem(f_shift,df)) > 1e-5
  error('The frequency shift is not a multiple of the channel spacing');
end 
%
ch_shift = round( f_shift / df );


%=== Number of input and output frequenies
%
nin   = length( f_sensor );
nout  = nin - ch_shift;


%=== 
if exist('rweights') & rweights
  ws = [-1;1];
else
  ws = [1;-1];
end


%=== Create transfer matrix for the frequency shifting, Hfs, and new F_SENSOR
%
rows = repmat( 1:nout, 2, 1);
rows = rows(:);
%
cols = [1:nout;(1+ch_shift):nin];
%
f_sensor = vec2col( f_sensor );
f_sensor = mean( f_sensor(cols) )';
%
cols = cols(:);
%
wgts = repmat( ws, nout, 1 );
%
Hfs = sparse( rows, cols, wgts, nout, nin );


%=== Include Hfs in H
%
H = h_x_h(Hfs,H);


%=== Create new F_Y and ZA_Y
%
[f_y,za_y] = h_fix_ys(f_sensor,za_sensor);
