%------------------------------------------------------------------------
% NAME:    h_x_h
%
%          Multiplicates two H matrices, H = Hpart*H, and checks if it is
%          more efficient to make H a full matrix instead of sparse.
%
% FORMAT:  H = h_x_h(Hpart,H)
%
% RETURN:  H        the new H matrix
% IN:      Hpart    the H matrix to be included in H
%          H        old H matrix
%------------------------------------------------------------------------

% HISTORY: 25.08.00  Created by Patrick Eriksson. 


function H = h_x_h(Hpart,H)


if size(H,1)==1 & size(H,2)==1
  if H ~= 1
    error('When H is a scalar, only the value 1 is allowed')
  end
  H = Hpart;

elseif size(Hpart,1)==1 & size(Hpart,2)==1
  if Hpart ~= 1
    error('When Hpart is a scalar, only the value 1 is allowed')
  end
  return

else
  if size(Hpart,2) ~= size(H,1)
    error('The two H matrices have incompatible sizes');
  end
  H = Hpart*H; 

end

clear Hpart


%=== Make H full if more than 2/3 of the values are used
if issparse(H) & (length(find(H)) > size(H,1)*size(H,2)*0.67)
  H = full(H);
end
