%------------------------------------------------------------------------
% NAME:    create_binfile
%
%          Generates a binning file. The binning is based on the derivatives 
%          in the spectra, a low intensity difference between neighbouring
%          channels results in a high binning and vice versa.
%          NOUT can be a scalar or a vector. If NOUT is a scalar, the number
%          of data points after binning will be the same for all zenith 
%          angles. If NOUT is a vector it must have the same length as ZA
%          and gives then the number of bins for each zenith angle.
%          The numbers in NOUT can be 0.
%          Channels with an optical thickness > TAU_LIM are removed and if
%          the number of remianing channels < NOUT that data from that
%          zenith angle is totally removed. 
%
% FORMAT:  create_binfile(f,za,tb,tau,nout,tau_lim,filename)
%
% RETURN:  -
% IN:      f            frequency grid
%          za           zenith angle grid
%          tb           intensity vector
%          tau          vector with optical thicknesses
%          nout         number of output values / zenith angles
%          tau_lim      limit on the optical thickness
%          filename     full name on output file
%------------------------------------------------------------------------

% HISTORY:     1997  A first version created for Skuld by Patrick Eriksson. 
%          00.08.30  Modified and improved for ARTS by Patrick Eriksson. 


function create_binfile(f,za,tb,tau,nout,tau_lim,filename)


f    = vec2col(f);
za   = vec2col(za);
nout = vec2col(nout);


%=== Main sizes
nf    = length(f);
nza   = length(za);


%=== Check input 
if length(nout) == 1
  nout = nout*ones(nza,1);
else
  if length(nout) ~= nza
     error('When NOUT is a vector, ZA and NOUT must have the same length');
  end
end
if length(tb) ~= nf*nza
   error('Lengths of F, ZA and TB do not match');
end
if length(tau) ~= nf*nza
   error('Lengths of F, ZA and TAU do not match');
end


%=== Convert TB and TAU to matrices
Tb    = reshape(tb,nf,nza);
Tau   = reshape(tau,nf,nza);


%=== Calculate a derivative measure
dTb       = [abs(diff(Tb))./(diff(f)*ones(1,nza));zeros(1,nza)];
dTb(nf,:) = dTb(nf-1,:);
dTb(2:nf-1,:) = (dTb(1:nf-2,:)+dTb(2:nf-1,:))/2;


%== Set the derivative to INF for optically thick channels
dTb(find(Tau>tau_lim)) = Inf;


%== Find number of output angles and set up output cell
nza_out = sum( (nout'>0) & (sum(~isinf(dTb))>nout') );
A       = cell(nza_out,1);


%=== Make heading text
heading = str2mat(...
   'BINNING FILE',...
   'This binning pattern has been created by the file create_binfile',...
   sprintf('Limit on optical thickness: %.2f',tau_lim),...
   sprintf('Number of input angles: %d',nza),...
   sprintf('Number of output angles: %d',nza_out),...
   sprintf('Number of input frequencies: %d',nf));
if length(nout) == 1
  heading = str2mat(heading,sprintf('Number of output frequencies: %d',nout));
else
  heading = str2mat(heading,...
          sprintf('Number of output frequencies: %d-%d',min(nout),max(nout)));
end


%=== loop spectra
for i =1:nza
  
  %= Find index of Inf
  i_inf = find(isinf(dTb(:,i)));
  n_inf = length(i_inf);

  if (nout(i)<1) | (n_inf> (nf-nout(i)))
     A{i} = 0;


  else

    start  = (1:nf)';
    chs    = ones(nf,1);
    w      = dTb(:,i);
    l      = nf;

    while l > (nout(i)+n_inf)

      ws = sum([w(1:l-1)';w(2:l)'])/2;

      [u,i_min] = min(ws);
      i_min     = i_min(1);

      start = start([1:i_min,i_min+2:l]);

      chs(i_min) = chs(i_min) + chs(i_min+1);
      chs = chs([1:i_min,i_min+2:l]);

      w(i_min) = w(i_min) + w(i_min+1);
      w = w([1:i_min,i_min+2:l]);
      
      l = l - 1;

    end

    ind = find(~isinf(w));
    start = start(ind);
    chs   = chs(ind);

    A{i} = [start start+chs-1];


  end % else

end %for i


write_datafile(filename,A,heading,0)
