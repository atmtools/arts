%------------------------------------------------------------------------
% NAME:    hEigs
%
%          Includes into H and Hd data reduction based on decomposition in 
%          eigenvectors.
%          The eigenvectors shall be stored as a H matrix with hSave,
%          where each column is an eigenvector.
%          The number of eigenvectors to use is controlled by the
%          parameter N. The saved data must of course have more than
%          N columns.
%          Note that the reduction matrix shall be saved as H. When
%          reading the data (by hLoad), Hd, f_sensor and za_sensor are
%          ignored. The first N values of F_Y and ZA_Y are returned.
%
% FORMAT:  [H,f_y,za_y,Hd] = hEigs(H,Hd,basename,n)
%
% RETURN:  H           total H matrix after data reduction
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          Hd          data reduction H matrix after data reduction
% IN:      H           total H matrix before data reduction
%          Hd          Hd matrix before data reduction
%          basename    basename for H file 
%          n           number of eigenvectors to use
%------------------------------------------------------------------------

% HISTORY: 00.08.25  Created by Patrick Eriksson. 


function [H,f_y,za_y,Hd] = hEigs(H,Hd,basename,n)


%=== Load the H matrix and associated grids
[Hnew,f_y,za_y] = hLoad(basename);


%=== Check data
if size(Hnew,2) > n
  error(sprintf('The saved data have not sufficent columns (%d, not %d)',size(Hnew,2),n));
end


%=== Perform multiplication
Hnew = Hnew(:,1:n)';
H    = h_x_h(Hnew,H);
Hd   = h_x_h(Hnew,Hd);


%=== Cut F_Y and ZA_Y
f_y  = f_y(1:n);
za_y = za_y(1:n);
