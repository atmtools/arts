%------------------------------------------------------------------------
% NAME:    rand_test
%
%          Plots the histogram and the autocorrelation function for the
%          data in a random vector. The main use of the function should
%          be to check if the random number generator fulfils the expected
%          criteria of PDF and uncorrelated data.
%
%          The expected PDF is plotted. As uncorrelated data is expected
%          the wanted autocorrleation function is a delta function at zero
%          lag.
%
%          The input data must match the expected PDFs. The following
%          distribution types are handled:
%
%          Uniform: TYPE='u', uniform PDF between 0 and 1
%                   To produce such a vector in ARTS:
%                      VectorRandUniform(z_abs){
%                                                low  = 0
%                                                high = 1
%                                                n    = 1000000 }
%                      VectorWriteBinary(z_abs){""}
%
%          Normal:  TYPE='n', normal PDF with mean 0 and std. dev. 1
%
% FORMAT:  rand_test(z,type,nbins,maxlag)
%
% RETURN:  -
% IN:      z        random vector of any length
%          type     PDF type:
%                     'u' for uniform
%                     'n' or 'g' for normal
%          nbins    number of bins for histogram
%          maxlag   maximum lag for autocorrelation function
%------------------------------------------------------------------------

% HISTORY: 00.11.27  Created by Patrick Eriksson. 


function rand_test(z,type,nbins,maxlag)

type = lower(type);


%=== Make histogram with relative frequencies
figure(1),clf
[n,x] = hist(z,nbins);
n = n/length(z);
bar(x,n,1);
xlabel('Normalised value');
ylabel('Relative frequency')
hold on
ax = axis;


%=== Plot expected distribution
if strcmp(type(1),'u')
  h = plot([0 1],[1 1]*(x(2)-x(1)),'r-');
elseif strcmp(type(1),'g') | strcmp(type(1),'n') 
  h = plot(x,normpdf(x,0,1)*(x(2)-x(1)),'r-');
else 
  error('Unknown PDF');
end
set(h,'LineWidth',2)
title('The red line is the expected PDF')



%=== Calculate and plot autocorrelation
figure(2),clf
%= Remove expected mean for uniform PDF
if strcmp(type(1),'u')
  z = z - 0.5;
end
[c,lags]=xcorr(z,maxlag,'coeff');   
plot(lags,c)
xlabel('Lag')
ylabel('Autocorrelation')
title('A delta function at 0 is expected')
