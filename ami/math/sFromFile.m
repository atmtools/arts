%------------------------------------------------------------------------
% NAME:    sFromFile
%         
%          Derives a covariance matrix by using the info given in file
%          as detailed below.
%
% FORMAT:         S = sFromFile(file,vector)
%
% RETURN:  S      Covariance matrix
%
% IN:      file   File where to get the s info. It should follow
%                 arts format as detailed below, and only ASCII files
%                 are allowed.
%
%                 The abscissa of the covariance matrix is determined
%                 by the vector. The covariance matrix can be 
%                 constructed as a sum of an arbitrary number of covariance
%                 matrices.
%
%                 The file starts with a matrix defining some parameters.
%                 The first row of the first matrix gives the
%                 correlation function flag for each covariance matrix,
%                 the possibilities are:
%                   0  no correlation, diagonal matrix. The corrlength
%                      and cutoff are of no importance here (see below)
%                   1  linearly decreasing to 0 (tenth function)
%                   2  exponential
%                   3  gaussian            
%                 The second row of the first matrix gives the correlation
%                 cut-off for each covariance matrix. Correlations below
%                 this value are set to 0. This variable can be used to make
%                 S more sparse.    
%                 Notice that the number of columns of the first matrix
%                 must equal the length of the matrix array - 1.
%
%                 Next are the covariance matrix/ces, each one is defined
%                 by a three column matrix, where column 1 is the abscissa
%                 (matching the given vector), column 2 standard deviations
%                 and column 3 correlation lengths.
%                 Linear interpolation is used to get values at intermediate
%                 points.
% 
%                 An example, where a diagonal matrix and a gaussian matrix,
%                 both having a standard deviation of 1, are summed, 
%                 the gaussian correlation length is 2 and k_grid has values
%                 between 0 and 10:
%                 3
%                 2 2
%                 0 2
%                 0 0
%                 2 3
%                 0   1  0
%                 10  1  0
%                 2 3
%                 0   1  2
%                 10  1  2             
% 
%         vector  Vector determining the abscisa of the covariance matrix.
%
%--------------------------------------------------------------------------
%
% HISTORY: 01.08.14  Created by Carlos Jimenez/Patrick Eriksson


function S = sFromFile( file, vector )


out(3,'Setting up a covariance matrix');


%=== Reading file
%
I = read_datafile(file,'AOMATRIX');


%=== Some checks of sizes
%
n = size(I,1);
%
if n < 1 
  error('The file must contain > 1 matrix.');
%
elseif size(I{1},1)~= 2 
  error('The first matrix in the file must have 2 rows.');
%
elseif size(I{1},2) ~= n-1 
  error( ['The number of columns of the first matrix must equal the',...
            'number of matrices - 1.'] );
end


%=== Loop the different covariance matrices
%
out(3,['Summing ',num2str(n-1),' matrices.']);
%  
for i=2:n

    % --  Check if the corrfun flag is an integer
    if (I{1}(1,i-1)-floor(I{1}(1,i-1))) ~= 0 
      error('The first row of matrix 1 shall only contain integers.');
    end

    % --  Move definition values to vectors
    kp      = I{i}(:,1);
    sdev    = I{i}(:,2);
    clength = I{i}(:,3);

    % -- Create covariance matrix
    if i == 2
      S = setup_covmatrix( vector, I{1}(1,i-1), I{1}(2,i-1), kp, sdev, clength );     
    else
      S = S + setup_covmatrix( vector, I{1}(1,i-1), I{1}(2,i-1), kp, sdev, clength );
    end

end

if isempty(S)
   error('Covariance matrix empty, check your grid vector')
end


%=== Make S sparse if less than 2/3 of the values are used
if length(find(S)) < size(S,1)*size(S,2)*0.67

     S = sparse(S);
end





%----------------------------------------------------------------------
%
%                              SUB-FUNCTIONS
%
%----------------------------------------------------------------------

%-------------------------------------------- ----------------------------
function  s = setup_covmatrix( kg, corrfun, cutoff, kp, sdev, clength ) 
%------------------------------------------------------------------------

%=== Checking
%
n = length(kg);
%
if length(sdev) ~= length(clength)
  error(['The standard deviation and correlation length vectors',...
         'must have the same length']);
%
elseif min(kg)<min(kp) | max(kg)>max(kp)
  error(['The data defining the covariance do not cover all',...
         'retrieval/error points.']);
%
elseif cutoff >= 1
  error('The correlation cut-off must be < 1.');
end




%=== Interpolate to get standard deviation and correlation length at
%    each point of k_grid
%
sd  = interp1( kp, sdev, kg );
cl  = interp1( kp, clength, kg );


%=== Diagonal matrix
%
if corrfun == 0
  %
  out(3,'Creating a diagonal covariance matrix');
  s = diag( sd.*sd );


%===  Linearly decreasing (tenth function)
%  
elseif corrfun == 1
  %  
  out(3,'Creating a simple covariance matrix');
  out(3, '  Tenth function correlation'); 
  out(3,['  Correlation cut-off : ',num2str(cutoff)]);
  out(3,['  Size                : ',num2str(n)]);
  s = zeros(n,n);
  for row=1:n
    for col=row:n
      c = 1 - 2 * (1-exp(-1)) * abs( (kg(row)-kg(col))/(cl(row)+cl(col)) );
      if c>0 & c>cutoff
        s(row,col) = c*sd(row)*sd(col); 
        s(col,row) = s(row,col);
      end
    end
  end


%===  Exponential
%  
elseif corrfun == 2
  %
  out(3,'Creating a simple covariance matrix');
  out(3, '  Exponential correlation');
  out(3,['  Correlation cut-off : ',num2str(cutoff)]);
  out(3,['  Size                : ',num2str(n)]);
  s = zeros(n,n);
  for row=1:n
    for col=row:n
      c = exp(-abs(kg(row)-kg(col))/((cl(row)+cl(col))/2));
      if c > cutoff
        s(row,col) = c*sd(row)*sd(col); 
        s(col,row) = s(row,col);
      end
    end
  end


%===  Gaussian
  
elseif corrfun == 3
  %
  out(3,'Creating a simple covariance matrix');
  out(3, '  Gaussian correlation');
  out(3,['  Correlation cut-off : ',num2str(cutoff)]);
  out(3,['  Size                : ',num2str(n)]);
  s = zeros(n,n);
  for row=1:n
    for col=row:n
%      c = exp( -1 * ( (kg(row)-kg(col))/(cl(row)+cl(col)) /2  )^ 2  );
      c = exp( -4 * ( (kg(row)-kg(col))/(cl(row)+cl(col)) )^2  );
      if c > cutoff
        s(row,col) = c*sd(row)*sd(col); 
        s(col,row) = s(row,col);
      end
    end
  end


%=== Unknown correlation function
  
else

  error(['Unknown correlation function flag ', num2str(corrfun)]);

end
