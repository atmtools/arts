%-----------------------------------------------------------------------------
% NAME:     test_Kt
%
%           A function to test if the temperature weighting functions are
%	    correct.
%
%	    The fields of Q shall be set in such way that spectra are
%           calculated (by qp_y) following the assumptions for the
%	    temperature weighting functions. Important is to match the
%           settings for hydrostatic equilibrium and refraction with the
%           settings for the weighting function calculation.
%
%	    A perfect match between dY and Kt shall not be expected but the
%           results shall be (very?) similar.
%
% FORMAT:   [Kt,dY] = test_Kt(Q)
%
% OUT:      Kt   Weighting function matrix for temperature.
%           dY   Disturbance calculations that should match Kt.
% IN:       Q    A Qpack structure.
%-----------------------------------------------------------------------------

% HISTORY: 2002-xx-xx  Created by Patrick Eriksson

% qp_H etc. should be run before calling the function

function [Kt,dY] = test_Kt(Q)


if Q.TEMPERATURE_DO ~= 3
  error('Temperature is not a retrieval variable.');
end


y = qp_y( Q );
[Dy,Kt,kx_names,kx_index] = qpcls( Q );


%= Extract Kt
ok = 0;
i  = 0;
while ~ok
  i = i + 1;
  if strncmp(kx_names(i,:),'Temperature:',12);
    ok = 1;
    ind = kx_index(i,1):kx_index(i,2);
  end
end
%
Kt = Kt(:,ind);


if nargout == 1
  return
end


%= Create a temporary directory
tmpdir    = temporary_directory( Q.TMP_AREA );


%= Copy species profiles and PTZ to tmp dir.
eval(['!cp ',Q.APRIORI_VMR,'* ',tmpdir])
eval(['!cp ',Q.APRIORI_PTZ,' ',fullfile(tmpdir,'ptz.aa')])


%= Set Q.APRIORI_VMR and Q.APRIORI_PTZ
[dummy,fname] = fileparts(Q.APRIORI_VMR);
Q.APRIORI_VMR = fullfile(tmpdir,fname);
Q.APRIORI_PTZ = fullfile(tmpdir,'ptz.aa');
%
ptz = read_datafile(Q.APRIORI_PTZ,'MATRIX');


%= Get temperature grid
grid = read_datafile(fullfile(Q.RETRIEVDEF_DIR,Q.TEMPERATURE_KGRID),'VECTOR');
n    = length(grid);


%= Loop and calculate spectra
dY = zeros( size(Kt) );


for i = 1:n

  dt = zeros(n,1);
  dt(i) = 1;

  t = ptz(:,2) + interpp(grid,dt,ptz(:,1));

  write_datafile(Q.APRIORI_PTZ,[ptz(:,1),t,ptz(:,3)],'MATRIX');

  dY(:,i) = qp_y( Q ) - y;

end


%=== DELETE THE TEMPORARY DIRECTORY
%
delete_tmp_dir( Q.TMP_AREA, tmpdir );
