%------------------------------------------------------------------------
% NAME:    hRedKx
%
%          Includes reduction based on Kx into Hd. The basic idea is to
%          do a Hotelling transformation of the spectral space. But as
%          the spectral variability can be decomposed as:
%
%                 Sy = KxSxKx' + KbSbKb' + Se
%
%          and we are interested in capture the variability of the 
%          parameters to be retrieved (x), the transformation is 
%          calculated only from that variability:  
%
%                 KxSxKx' = E V E'
%
%          where E are the eigenvectors and V the eigenvalues.
%
%          The Hotteling transformation is derived as:      
%          
%          1. Notice that
%            KxSxKx' = Kx sqrt(Sx) sqrt(Sx)' Kx' = Kx sqrt(Sx) (Kx sqrt(Sx))'
%          where sqrt(Sx) is one of the possible square roots of Sx.
%          
%          2. Then a SVD
%                    Kx sqrt(Sx) = U M V'   
%          where U and V are left and right singular vectors.
%
%          3. The eigenvectors are then U as from (1):
%                    E V E' = U M V' (U M V')' = U M M U'
%
%          Different assumptions are implemented for sqrt(Sx):
%
%          1.Q.KRED_DEPTH = 0 
%                    sqrt(Sx) = I
%          i.e. Sx is assumed to be the identity matrix
%  
%          2.Q.KRED_DEPTH = 1 
%                    sqrt(Sx) = sqrt(diag(Sx))
%          i.e. Sx is assumed to be diagonal
% 
%          1.Q.KRED_DEPTH = 2 
%                    sqrt(Sx) = chol(Sx)
%          i.e. full Sx is used and the sqrt is done by a choleski
%          decomposition.
% 
%          More info in:
%              Eriksson, P., Jimenez, C., Búhler, S. and Murtagh, D. "A 
%              Hotteling Transformation approach for rapid inversion of   
%              atmospheric spectra", J. Quant Rad. Spectroscopy, pp    
%              2001.
%
%
%
% FORMAT:  [H,f_y,za_y,Hd] = hRedKx(Q,Hd,f_y,za_y,Hd)
%
% RETURN:  Q           Q structure
%          H           total H matrix after binning
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          Hd          data reduction H matrix after binning
% IN:      H           total H matrix before reduction
%          f_y         input frequency vector 
%          za_y        input zenith angle vector
%          Hd          data reduction H matrix before reduction
%          
%
% HISTORY: 01.07.18  Created by Carlos Jimenez/ Patrick Eriksson 
%------------------------------------------------------------------------


function [H,f_y,za_y,Hd] = hRedKx(Q,H,f_y,za_y,Hd)


out(1,1);
out(1,'Setting up reduction matrix based on Kx (H)');

%=== Create cfile

QE.DO_NONLIN  = 0;
QE.RECALC_ABS = 0;
template      = which( 'oem_iter1.tmplt' );
tmpdir    = temporary_directory( Q.TMP_AREA );
[cfile,basename] = qtool( Q, tmpdir, template, QE );

%=== Run ARTS and load results
  
qp_arts( Q, cfile ); 
Kx = read_artsvar( basename, 'kx' );

Kx = h_x_h(H,Kx);

%=== Scale KX by elements of SX

if     Q.KRED_DEPTH == 2
  
  load( [Q.OUT,'.sx'], '-mat' );
  np = size( Sx, 1 );
  Kx = h_x_h(Kx,chol(Sx));

elseif Q.KRED_DEPTH == 1

  load( [Q.OUT,'.dxinv'], '-mat' );
  np = size( Dxinv, 1 );
  Kx = h_x_h(Kx,sqrt(inv(Dxinv)));

end

%=== Calculate the reduction transfer matrix
[U,S,V]    = svd( Kx, 0 );
clear S V Kx
if Q.KRED_N > size(U,2)
    error(sprintf( ...
   '%d number of eigenvectors cannot be provided.\nOnly %d are possible.',...
    						 Q.KRED_N,size(U,2)));
end
U    = U( :, 1:Q.KRED_N )';
H    = h_x_h( U, H );
Hd   = h_x_h( U, Hd );
clear U
f_y  = 1:Q.KRED_N;
za_y = ones( Q.KRED_N, 1 );  

%=== Delete the temporary directory
delete_tmp_dir( Q.TMP_AREA, tmpdir );

out(1,-1);
