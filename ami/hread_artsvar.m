%------------------------------------------------------------------------
% NAME:    hread_artsvar
%
%          Reads data from a ARTS output file and applies a H matrix.
%
%          The data is read from the file
%             basename.varname.am
%
%          See further READ_ARTSVAR
%
% FORMAT:  x = hread_arts(H,basename,varname)
%
% RETURN:  x           the data
% IN:      H           the H matrix to apply
%          basename    the ARTS basename 
%          varname     variable name
%------------------------------------------------------------------------

% HISTORY: 17.04.00  Created by Patrick Eriksson. 

function x = hread_artsvar(H,basename,varname)


%=== Read the data by using READ_ARTSVAR
x    = read_artsvar(basename,varname);


%=== Apply the H matrix if not a scalar
if (size(H,1)==1) & (size(H,2)==1)

  if H ~= 1
    error('When H is a scalar, only the value 1 is allowed')
  end

else

  %=== Check sizes
  if size(H,2) ~= size(x,1)
    error(sprintf('Sizes of H and the data do not match (%d and %d, respectively)',size(H,2),size(x,1)));
  end

  x    = H*x;

end
