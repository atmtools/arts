%------------------------------------------------------------------------
% NAME:    xmlLoad
%
%          Loads a xml file.
%
% FORMAT:  result = xmlLoad(filename)
%
% RETURN:  result      Data read from file
% IN:      filename    XML filename
%------------------------------------------------------------------------

% HISTORY: 2002-09-25  Created by Oliver Lemke.

function result = xmlLoad(filename)

fid = fopen (filename,'r');

if fid == -1
  error (sprintf ('Cannot open file %s', filename));
end

%=== Validate XML file header
s = fscanf (fid, '%s', 1);
limit = 10;
while limit & s(size(s)) ~= '>'
  s = strcat (s, fscanf (fid, '%s', 1));
end

if ~strcmp (s, '<?xmlversion="1.0"?>')
  error ('Invalid xml header');
end

%=== Parsing tags
while ~feof (fid)
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
      
      attrlist = xmlReadAttributes (fid);
      
      switch tag
       case 'arts'
       otherwise
	func = str2func (strcat ('xmlRead', tag));
	result = feval (func, fid, attrlist);
	disp (result);
      end
    else %=== or is it a closing tag
    end
  end
end

fclose (fid);

result = 0;
