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

%=== Parsing data tag
result = xmlReadTag(fid);

fclose (fid);

