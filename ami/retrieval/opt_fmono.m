%------------------------------------------------------------------------------
% NAME:     opt_fmono
%           
%           This is based on qpopt_fmono.m in Qpack/Setup.     
%
%           This function optimizes the monochromatic frequency grid based
%           on the polynomial order (F_ORDER).
%
%           The optimization is performed in this way:
%
%            1. Prepare spectra calculated on fine test grid.
%            2. The end points of the frequency range are included in
%               F_MONO.
%            3. The spectra are interpolated from the present F_MONO
%               to the fine test grid. Interpolation applied follows F_ORDER.
%            4. The maximum difference between interpolated and calculated
%               spectra (for the fine grid) is determined.
%            5. If the difference exceeds the given precision, the frequency
%               where the maximum difference was found is put into F_MONO
%            6. Points 3-5 are repeated until the wanted precision is obtained.
%
% FORMAT:   f_mono = opt_fmono( Y, ftest, prec, f_order 
%                                                  [, prec_unit, do_plot ])
%           
%
% OUT:      f_mono         Optimized monochromatic frequencies.
%
% IN:       Y              Spectra calculated on FTEST.
%                        
%                          [nf, nza, nsu] = size(Y)  
%
%                          nf  --- number of frequencies
%                          nza --- number of elements of za_pencil
%                          nsu --- number of set-up atmospheres  
%                                        
%           ftest          fine test frequency grid.  
%                         
%           prec           Required precision. The unit is either the 
%                          abosolute or the relative. 
%                          Default is the absolute unit. See below.
%
%           f_order        Interpolation order (1 or 3).
%
% OPTIONAL:
%
%           prec_unit      Flag for the unit of precision.
%                          Absolute unit  --- 1 (default)
%                          relative unit  --- 0 
%                          
%           do_plot        Flag to show example on position of found frequency
%                          grid and interpolation error. Default is 1 (meaning
%                          plotting). Note that [K] is used as a unit of Y.
%------------------------------------------------------------------------------

% HISTORY: 2001.11.20  Created by Patrick Eriksson.
%          2003.08.20  Modified by Sho Tsujimaru.  

function f_mono = opt_fmono( Y, ftest, prec, f_order, prec_unit, do_plot )

%== 
%
if ~exist('prec_unit'),  prec_unit = 1; end
if ~exist('do_plot'),    do_plot   = 1; end

%== Check F_ORDER
%
if f_order ~= 1 & f_order ~= 3
  error('F_ORDER must be 1 or 3.');
end

out(1,1);
%
out(1,'Optimising F_MONO with respect to interpolation error.');
%
if f_order == 1
  out(1,'Polynomial representation: linear');
else
  out(1,'Polynomial representation: cubic');
end
%
if prec_unit
  out(1,sprintf('Required precision (absolute) %.5f',prec));
else
  out(1,sprintf('Required precision (relative) %.5f',prec));
end 
%
if do_plot
  out(1,'Plotting will be performed.');
end
%
out(2,0);

%
f_mono = [];
%

% Size of Y: frequnecy, za_pencil, set-up atmosphere  
[nf, nza, nsu] = size( Y );

if nf ~= length(ftest)
  error('The length of FTEST is not consistent to that of Y.');
end

%
if f_order == 1
  ind  = 1;
else
  ind  = [1;round(nf/2)];
end
  
%
imax = nf;
%
emax = 1e99;
%

while emax > prec
  %
  ind = [ind;imax];
  ind = sort( ind );
  %
  emax = 0;
  %
  for j = 1:nsu
    if f_order == 1
      Y2 = interp1( ftest(ind), Y(ind,:,j), ftest, 'linear' );
    else
      Y2 = interp1( ftest(ind), Y(ind,:,j), ftest, 'cubic' );
    end
    % 
    if prec_unit
      [ej,ij] = max( max( abs(Y2-Y(:,:,j)), [], 2 ) ); 
    else 
      [ej,ij] = max( max( abs( (Y2-Y(:,:,j))./Y(:,:,j) ), [], 2 ) );
    end
    %
    if ej > emax
      emax = ej;
      imax = ij;
    end
  end
  %
  if emax > prec
    out(3,sprintf('Error=%.5f, adding %.6f GHz, nr=%d',...
                                        emax, ftest(imax)/1e9,length(ind)+1));
  else
      out(3,sprintf('Error=%.5f',emax));
  end
end

%= Show results if DO_PLOT
%
if do_plot
  %
i1 = ceil(nza/2); 
  %
  figure(1), clf
  plot( ftest/1e9, Y(:,i1,nsu), '-', ftest(ind)/1e9, Y(ind,i1,nsu), '.' );
  xlabel('Frequency [GHz]');
  ylabel('Tb [K]');
  title( sprintf( ...
        'Last test atmosphere, %dth zenith angle', i1));
  %
  figure(2), clf
  xlabel('Frequency [GHz]');
  % 
  if prec_unit 
    plot( ftest/1e9, (Y2(:,i1) - Y(:,i1,nsu))*1e3, '-' );
    ylabel('Interpolation error [mK]');
  else
    plot( ftest/1e9, (Y2(:,i1) - Y(:,i1,nsu))./Y(:,i1,nsu)*1e2, '-' );
    ylabel('Interpolation error [%]');
  end 
  % 
  title( sprintf( ...
        'Last test atmosphere, %dth zenith angle', i1));
  %
end

f_mono = [ f_mono; ftest(ind) ];

%=== Sort and remove duplicates of frequencies
%
f_mono = unique( f_mono );

out(2,0);
out(1, sprintf('No.of F_MONO before optimization = %d', nf));
out(1, sprintf('No.of F_MONO after optimization = %d', length(f_mono)));

out(1,-1);
