%------------------------------------------------------------------------
% NAME     filedate
%
%          Returns the modification date for a file, either as a date
%          string or as a serial date number.
%          To check if a file is newer than another
%          FILEDATE(filename1,1) < FILEDATE(filename2,1) ...
%
% RETURN   date      the date
% IN       filename  the filename
% OPTIONAL asnumber  ASNUMBER=0, date is a string (default)
%                    ASNUMBER=1, date is a number
%------------------------------------------------------------------------

% HISTORY: 99.10.13   Created by Patrick Eriksson.


function date = filedate(filename,asnumber)


if ~exist('asnumber')
  asnumber = 0;
end


d    = dir(filename);
date = d.date;


if asnumber
  date = datenum(date);
end
