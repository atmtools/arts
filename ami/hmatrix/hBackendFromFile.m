%------------------------------------------------------------------------
% NAME:    hBackendFromFile
%
%          Includes into H a backend where the channel response is defined 
%          in a file. Channels can (1) have the same response, or (2) have
%          individual responses. The channel file shall have ARTS format. For 
%          case common repsonse, it shall be a 2 column matrix, where column 1
%          is (relative) frequencies and column 2 is response values. For 
%          individual resposnse channel, it shall be an array of matrices, 
%          where matrix {i} shall be a 2 column matrix giving the response 
%          of channel {i}, also column 1 frequencies and column 2 response 
%          values  The response of the channels is normalised and the
%          response values do not need to be normalised.
%
% FORMAT:  [H,f_y,za_y,f_sensor] = hBackendFromFile(H,f_sensor,za_sensor,
%                                                f_centre,filename,o_ch,o_y)
%
% RETURN:  H           H matrix after backend
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          f_sensor    new frequency grid
% IN:      H           H matrix before the backend
%          f_sensor    input frequency grid
%          za_sensor   zenith angles
%          f_centre    centre frequencies of the channels
%          filename    name on file with channel response specification
%          o_ch        linear (=1) or cubic (=3) treatment of the channel
%                      response
%          o_y         linear (=1) or cubic (=3) treatment of spectra
%------------------------------------------------------------------------

% HISTORY: 00.08.25  Created by Patrick Eriksson. 
%          00.11.16  Included linear/cubic flags (PE) 
%          05.05.26  Modified to have the possibility of
%                    different channel response (CJ ) 

function [H,f_y,za_y,f_sensor] = ...
             hBackendFromFile(H,f_sensor,za_sensor,f_centre,filename,o_ch,o_y)


%=== Read the file defining the channel response
%    first, it looks for a backend having same response
%    for all channels, if it fails it looks for a
%    backend having individual responses 

try

  A = read_datafile( filename, 'MATRIX' );
  Af = A(:,1);
  Ab = A(:,2);

catch

  A = read_datafile( filename, 'AOMATRIX' );
  for j = 1:length( A )
    Af{j} = A{j}(:,1);
    Ab{j} = A{j}(:,2);
  end

end



%=== Get H for the backend
[Hback,f_sensor] = h_backend(f_sensor,f_centre,za_sensor,Af,Ab,o_ch,o_y);


%=== Include Hback in H
H = h_x_h( Hback, H );


%=== Create new F_Y and ZA_Y
[f_y,za_y] = h_fix_ys( f_sensor, za_sensor );
