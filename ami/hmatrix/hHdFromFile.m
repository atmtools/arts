%------------------------------------------------------------------------
% NAME:    hHdFromFile
%
%          Loads a data reduction H matrix (Hd) and includes this in the total
%          transfer and data reduction matrices (H and Hd, respectively).
%          The Hd matrix shall be saved with hSave (as Hd, the H matrix is
%          ignored, as F_SENSOR and ZA_SENSOR). The saved F_Y and ZA_Y are
%          returned.
%
% FORMAT:  [H,f_y,za_y,Hd] = hHdFromFile(H,Hd,basename)
%
% RETURN:  H           new total H matrix
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          Hd          new data reduction H matrix
% IN:      H           total H matrix before data reduction
%          Hd          old data reduction H matrix
%          basename    basename for H file 
%------------------------------------------------------------------------

% HISTORY: 25.08.00  Created by Patrick Eriksson. 


function [H,f_y,za_y,Hd] = hHdFromFile(H,Hd,basename)


%=== Load the Hd matrix and associated grids
[U,f_y,za_y,Hnew] = hLoad(basename);


%=== Perform multiplication
H    = h_x_h(Hnew,H);
Hd   = h_x_h(Hnew,Hd);



