%------------------------------------------------------------------------
% NAME:    hRedLimb
%
%          Includes reduction based on Limb into Hd. By Limb it means
%          a special reduction based on combining the reduction based
%          on Kx with a prior optimization of the spectral input 
%          from a limb observation regarding each specific point of a 
%          retrieval grid. This reduction is useful for a retrieval
%          method where each element of the retrieved species state 
%          vector is obtained independently, so the spectral input
%          in the limb observation can be optimized for each element.
%          So far we use it with the net retrievals from the package
%          Npack.
%
%          The basic idea is to decide which elements for the spectral
%          input are left for each retrieval point, and then to proceed
%          with a normal Hotelling transformation of the new spectral
%          space. In practice is done by doing the Hot transf in a Kx
%          padded with zeros to leave only the rows (part of the scan)
%          and columns (the retrieval tag) corresponding to each
%          retrieval tag and point.
%
%          Notice that now for each point of the retrieval grid there
%          will be a reduction matrix, also saved in Hd as described
%          below. The same applies to the subsequent H. 
%
%          The parameters for this reduction are given as usual as
%          part of Q, see README.    
%          
%          For details of the Hotelling transformation see hRedKx.
%
%          NOTE: Warning !! Under test !!
%
% FORMAT:  [H,f_y,za_y,Hd] = hRedLimb(Q,Hd,f_y,za_y,Hd)
%
% RETURN:  Q           Q structure
%          H           [] empty to signal in qp_H that H should not be
%                      saved again as it is saved here. [Q.OUT,'.h']
%                      will contain Hij matrices corresponding to
%                      the tag i and ret point j.
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          Hd          [] empty as in H, same saving strategy
%
% IN:      H           total H matrix before reduction
%          f_y         input frequency vector 
%          za_y        input zenith angle vector
%          Hd          data reduction H matrix before reduction
%          
% TEMPLATE:            oem_iter1.tmplt
%
% HISTORY: 01.08.14  Created by Carlos Jimenez
%------------------------------------------------------------------------


function [H,f_y,za_y,Hd] = hRedLimb(Q,H,f_y,za_y,Hd)


out(1,1);
out(1,'Setting up reduction matrix optimizing the limb scan (H)');

%=== Create cfile

out(2,'Running arts to get variables');
template      = which( 'oem_iter1.tmplt' );
tmpdir    = temporary_directory( Q.TMP_AREA );

QE.DO_NONLIN  = 0;
QE.RECALC_ABS = 0;
[cfile,basename] = qtool( Q, tmpdir, template, QE);

%=== Run ARTS and load results
  
qp_arts( Q, cfile );
%
f_mono     = read_artsvar( basename, 'f_mono' );
%
z_tan      = read_artsvar( basename, 'z_tan' );
p_abs      = read_artsvar( basename, 'p_abs' );
z_abs      = read_artsvar( basename, 'z_abs' );
za_pencil  = read_artsvar( basename, 'za_pencil' );
%
kx         = read_artsvar( basename, 'kx' );
kx_aux     = read_artsvar( basename, 'kx_aux' );
kx_names   = read_artsvar( basename, 'kx_names' );
kx_lengths = read_artsvar( basename, 'kx_lengths' );
%
nrq        = size( kx_names, 1 );
nf         = length(f_mono);
nza        = length(za_pencil);
kx_index   = [ [1;1+cumsum(kx_lengths(1:(nrq-1)))], cumsum(kx_lengths) ];

clear za_pencil kx_lengths 

%=== Setting new H and Hd
% 
%    To accomodate the new H and Hd - one for each
%    ret point and ret grid, matrices will be saved  
%    on [Q.OUT,'.h'] and [Q.OUT,'.hd'] as usual 
%    but with one matrix for tag and ret point
%    following
%      Hij 
%      Hdij   i=tag number  j=ret point


%=== Do data reduction from each point of each ret grid

out(2,'Doing data reduction');
   
h   = H;
hd  = Hd;
clear H Hd


% -- loading retrieval grid
p_gr    = find(Q.LRED_KGRIDS  =='"');
n_sp    = size(p_gr,2)/2;
for i=1:n_sp
  is = num2str(i);
  pa = num2str(p_gr(2*(i-1)+1)+1);
  pb = num2str(p_gr(2*(i-1)+2)-1);
  eval(['p_grid=Q.LRED_KGRIDS(',pa,':',pb,');'])
  f_grid=[Q.RETRIEVDEF_DIR,'/',p_grid];
  eval(['p_grid',is,'= read_datafile(''',f_grid,''',','''VECTOR''',');']);
end


for i = 1:nrq

  is     = num2str(i);      
  out(2,['Tag ',is]);
  ind    = kx_index(i,1):kx_index(i,2);
  eval(['p_grid = p_grid',is,';']);
  z_grid = interpp(p_abs,z_abs,p_grid);
  l_grid = length(z_grid);
 
  % -- leaving only eigenvectord for i species
  kxi = kx(:,ind);

  for j = 1:l_grid
 
      js    = num2str(j);
      out(2,['  Retrival altitude: ',num2str(z_grid(j)/1e3),' km']);  
      [h_n,hd_n] = limbKx(Q,h,hd,kxi,z_grid(j),z_tan,nf,nza);
      eval(['H',is,js,'=h_n;'])
      if i==1 & j==1
        eval(['save ',Q.OUT,'.h H',is,js]) 
      else
        eval(['save ',Q.OUT,'.h H',is,js,'   -append']) 
      end   
      eval(['clear H',is,js,' h_n;'])
      eval(['Hd',is,js,'=hd_n;'])
      if i==1 & j==1
        eval(['save ',Q.OUT,'.hd Hd',is,js]) 
      else
        eval(['save ',Q.OUT,'.hd Hd',is,js,' -append']) 
      end
      eval(['clear Hd',is,js,' hd_n;'])


  end

end

%===  Saving aux var nrq and l_grid
ind = [];
for i=1:nrq

  eval(['l = length(p_grid',is,');']);
  ind = [ind l];

end  
eval(['save ',Q.OUT,'.hd ind  -append'])   
eval(['save ',Q.OUT,'.h ind   -append'])   



%=== To signal to qp_H that H and Hd should not been
%    saved as they are saved here.

H    = [];
Hd   = [];    

%=== New f_y and za_y

f_y  = 1:Q.LRED_N;
za_y = ones( Q.LRED_N, 1 );  

%=== Delete the temporary directory

delete_tmp_dir( Q.TMP_AREA, tmpdir );
out(1,-1);

%----------------------------------------------------------------------
%
%                              SUB-FUNCTIONS
%
%----------------------------------------------------------------------

%
%-------------------------------------------- ----------------------------
  function [H,Hd] = limbKx(Q,H,Hd,Kx,z,ztan,nf,nza);
%------------------------------------------------------------------------


%=== Looking for z + HISTEP > z_tan > z - LOSTEP

 
hi       = Q.LRED_HISTEP;
lo       = Q.LRED_LOSTEP;
z_ind    = find( ( (z+hi)>ztan ) & ( ztan > (z-lo) ));     

% -- special cases

if isempty(z_ind)
  %
  %--  when z is lower or higher than any z_tan then
  %    kept as many scans as REDL_HISTEP or LOSTEP
  %    from the lowest or higher z_tan 

  if    z < ztan(length(ztan))
    %
    z_ind = find( (ztan(length(ztan))+hi)>ztan );
    %
  elseif z > ztan(1)
    %
    z_ind = find( ztan > (ztan(1)-lo));
    %
  end
  %
end


%=== Building the Kx corresponding to the optimized grid


% -- 1st leaving only the Kx corresponding to the nrq species 
%    (already done here)


% -- 2nd leaving only the Kx right channels
%    (padding with zeros the rows which channels are excluded)


nind     = length(z_ind);  
ind      = (1+ nf * (z_ind(1)-1)):(nf * z_ind(length(z_ind)));

Ka       = Kx(ind,:);
Kx       = zeros(size(Kx));
Kx(ind,:)= Ka;


%=== Apply sensor to Kx

Kx = h_x_h(H,Kx);


%=== Scale Kx by elements of Sx

if     Q.LRED_DEPTH == 2
  
  load( [Q.OUT,'.sx'], '-mat' );
  np = size( Sx, 1 );
  Kx = h_x_h(Kx,chol(Sx));

elseif Q.LRED_DEPTH == 1

  load( [Q.OUT,'.dxinv'], '-mat' );
  np = size( Dxinv, 1 );
  Kx = h_x_h(Kx,sqrt(inv(Dxinv)));

end


%=== Calculate the reduction transfer matrix
[U,S,V]    = svd( Kx, 0 );
clear S V Kx
if Q.LRED_N > size(U,2)
    error(sprintf( ...
   '%d number of eigenvectors cannot be provided.\nOnly %d are possible.',...
    						 Q.LRED_N,size(U,2)));
end


%=== saving results
U    = U( :, 1:Q.LRED_N )';
H    = sparse(H);
H    = h_x_h( U, H );
Hd   = sparse(Hd);
Hd   = h_x_h( U, Hd );
