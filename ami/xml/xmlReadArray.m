%------------------------------------------------------------------------
%
% NAME:    xmlReadArray
%
%          Reads Array from XML file.
%
% FORMAT:  result = xmlReadArray(fid, attrlist)
%
% RETURN:  result      Data read from file
% IN:      fid         File descriptor
% IN:      attrlist    List of attributes
%------------------------------------------------------------------------
%

% HISTORY: 2002-09-25  Created by Oliver Lemke.

function result = xmlReadArray(fid, attrlist)

ne = str2num (xmlGetAttrValue (attrlist, 'nelem'));

for e = 1:ne
  fprintf ('Reading array element');
  s = fscanf (fid, '%s', 1);
  
  if ~size (s)
    break
  end
  
  %=== Tag has to start with bracket
  if s(1) == '<'
    %=== Do we have an opening tag here?
    if s(2) ~= '/'
      l = size(s);
      tag = s(2:l(2));
      
      attrlist2 = xmlReadAttributes (fid);
    
      fprintf (tag);
      switch tag
       otherwise
	func = str2func (strcat ('xmlRead', tag));
	result{e} = feval (func, fid, attrlist2);
      end
    else %=== or is it a closing tag
    end
  end
end
