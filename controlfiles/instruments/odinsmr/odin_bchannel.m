% ODIN_BCHANNEL   Odin-SMR backend channel response in ARTS format
%
%     Returns data that match ARTS' *backend_channel_response*.
%
%    If stored to file:
%      xmlStore( filename, R, 'ArrayOfGriddedField1' );
%
% OUT   R        The response, as a ArrayOfGriddedField1.
% IN    df       Frequency seperation between channels.
% OPT   nohann   Set t true of no Hanning is appled. Default is false.

% 2012-11-17 Patrick Eriksson


function R = odin_bchannel(df,nohann)
%  
if nargin < 2
  nohann = false;
end
rqre_datatype( df, {@istensor0} );
rqre_in_range( df, 125e3, 1e6 );
rqre_datatype( nohann, {@isboolean} );


if nohann
  f = -22 : 0.5 : 22;
  r = sin(f) ./ f;
  r(isnan(r)) = 1;
  dd = fwhm( f, r )
  f  = (df/dd) * f;
else
  D = refdata;
  f = D(:,1) * df / 1e6;
  r = D(:,2);
end

R.name      = 'Backend channel response function';
R.gridnames = { 'Frequency' };
R.grids     = { f };
R.dataname  = 'Response';
R.data      = r;

% Convert to an array
%
R = { R };

return


function D = refdata

D = [
-3.937500000e+06 1.087367268e-03 
-3.718750000e+06 5.157498099e-03 
-3.500000000e+06 8.084014476e-03 
-3.281250000e+06 7.678047823e-03 
-3.062500000e+06 2.420031463e-03 
-2.843750000e+06 -7.445383788e-03 
-2.625000000e+06 -1.901840364e-02 
-2.406250000e+06 -2.642741522e-02 
-2.187500000e+06 -2.135781879e-02 
-1.968750000e+06 5.510310117e-03 
-1.750000000e+06 6.235955541e-02 
-1.531250000e+06 1.538415644e-01 
-1.312500000e+06 2.790393242e-01 
-1.093750000e+06 4.303879152e-01 
-8.750000000e+05 5.939781276e-01 
-6.562500000e+05 7.513495512e-01 
-4.375000000e+05 8.825018651e-01 
-2.187500000e+05 9.695182402e-01 
0.000000000e+00 1.000000000e+00 
2.187500000e+05 9.695182402e-01 
4.375000000e+05 8.825018651e-01 
6.562500000e+05 7.513495512e-01 
8.750000000e+05 5.939781276e-01 
1.093750000e+06 4.303879152e-01 
1.312500000e+06 2.790393242e-01 
1.531250000e+06 1.538415644e-01 
1.750000000e+06 6.235955541e-02 
1.968750000e+06 5.510310117e-03 
2.187500000e+06 -2.135781879e-02 
2.406250000e+06 -2.642741522e-02 
2.625000000e+06 -1.901840364e-02 
2.843750000e+06 -7.445383788e-03 
3.062500000e+06 2.420031463e-03 
3.281250000e+06 7.678047823e-03 
3.500000000e+06 8.084014476e-03 
3.718750000e+06 5.157498099e-03 
3.937500000e+06 1.087367268e-03 
];
return