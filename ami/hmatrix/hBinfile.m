%------------------------------------------------------------------------
% NAME:    hBinfile
%
%          Includes binning of spectromer channels into Hd. That is,
%          averages of neighbouring channels are used.
%          The binning is specified by a ARTS file with the same number of
%          matrices as viewing angles before binning. The matrices have
%          two columns where the first and second column gives the first 
%          and last channel, respectively, to join. 
%          If the first element of the matrix (element (1,1)) is zero,
%          the data from that viewing angle is totally neglected. 
%          The number of values after binning can vary between the 
%          zenith anngles.
%          An example:
%
%            # Number of matrices
%            3
%            #
%            # Viewing angle 1
%            2 2
%            #
%            1 2
%            3 3
%            #
%            # Viewing angle 2
%            1 2
%            #
%            1 3
%            # Viewing angle 3
%            1 1
%            #
%            0
%
%          For angle 1, channel 1 and 2 are combined while channel 3 is 
%          untouched. For angle 2, all three channels are combined.
%          Angle 3 is removed.
%          Note that the channels can be picked out in any order. In
%          addition, a channel can be used for several averages (but
%          make little sense).           
%
% FORMAT:  [H,f_y,za_y,Hd] = hBinFile(Hd,f,za,filename)
%
% RETURN:  H           total H matrix after binning
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          Hd          data reduction H matrix after binning
% IN:      H           total H matrix before binning
%          f_y         input frequency vector 
%          za_y        input zenith angle vector
%          Hd          data reduction H matrix before binning
%          filename    name on file specifying the binning pattern
%------------------------------------------------------------------------

% HISTORY: 00.08.22  Created by Patrick Eriksson. 


function [H,f_y,za_y,Hd] = hBinFile(H,f_y,za_y,Hd,filename)

disp('WARNING, this function does so far not consider the channel widths');


%=== Check F_Y and ZA_Y
ny  = length(f_y);     % total length
if ny ~= length(za_y)
  error('F_Y and ZA_Y must have the same length');
end


%=== Load binning file
bdata = read_datafile(filename);
% Conversion needed for Matlab
if ~iscell(bdata)
  temp = bdata;
  bdata = cell(1);
  bdata{1} = temp;
  clear temp
end
nb    = length(bdata);


%=== Some checks if binning pattern and get number of input and output values
nza  = 0;          % Number of output angles
nout = 0;          % Number of output values
nin  = 0;          % Number of input values used
for i = 1:nb
  if bdata{i}(1,1) > 0
    if size(bdata{i},2) ~= 2
      error('The binning patterns must be given as 2 column matrices');
    end
    if min(bdata{i}(:,2)-bdata{i}(:,1)) < 0
      error(sprintf('The value in column 2 cannot be smaller than the one in column 1 (angle %d)',i));
    end
    if min(min(bdata{i}))<1
      error(sprintf('Binning value(s) < 1 for angle %d',i));
    end
    nza  = nza + 1;
    n    = size(bdata{i},1);
    nout = nout + n; 
    nin  = nin + sum(bdata{i}(:,2)-bdata{i}(:,1)) + n;
  end
end


%=== Set up F_Y and ZA_Y and some vectors defining the local Hd
f_y2   = zeros(nout,1);          % new frequency vector
za_y2  = zeros(nout,1);          % new zenith angle vector
row    = zeros(nin,1);           % row index in Hd
col    = zeros(1,nin);           % col index in Hd
w      = zeros(nin,1);           % the values in Hd
iout   = 0;                      % output index
iin    = 0;                      % input index
iza    = 0;                      % index for zenith angles
ilim   = [0 0];                  % first and last index/col for iza
for i = 1:nb

  iza = iza + 1;
  ilim(1) = ilim(2)+1;
  if ilim(1) > ny
    error('Your binning file has too many zenith angles');
  end
  ilim(2) = ilim(1);
  while (ilim(2)<ny) & (za_y(ilim(2)+1)==za_y(ilim(1)))
    ilim(2) = ilim(2) + 1;
  end

  if bdata{i}(1,1) > 0
    for j = 1:size(bdata{i},1)

      i1            = bdata{i}(j,1);
      i2            = bdata{i}(j,2);
      iout          = iout + 1;
      wv            = 1/(i2-i1+1);
      
      if (ilim(1)+i2-1) > ilim(2)
        error(sprintf('You have selected value %d, but zenith angle %d has only %d values',i2,i,ilim(2)-ilim(1)+1));
      end

      for k = i1:i2
        iin         = iin + 1;
        row(iin)    = iout;
        col(iin)    = ilim(1)+k-1;
        w(iin)      = wv;
      end      

      i1            = ilim(1)+i1-1;
      i2            = ilim(1)+i2-1;
      f_y2(iout)    = (f_y(i1)+f_y(i2))/2;
      za_y2(iout)   = za_y(i1);
    end
  end
end


%=== Final check of number of zenith angles
if ilim(2) ~= ny
  error('Your binning file has too few zenith angles');
end


%=== Include the binning in H and Hd
Hb = sparse(row,col,w);
H  = h_x_h(Hb,H);
Hd = h_x_h(Hb,Hd);


%=== Copy F_Y2 to F_Y etc
f_y  = f_y2;
za_y = za_y2;
