%------------------------------------------------------------------------
% NAME:     vtag2arts
%
%           Converts a Verdandi tag number to the corresponding species
%           name according to the ARTS convention. A Verdandi tag number
%           is the molecule and isotope numbers joined, e.g. 31 for the
%           most common ozone isotope.
%
%           The Verdandi tag numbers are an extension of the HITRAN 
%           ones, so this function should also work for HITRAN tag
%           numbers.
%
%           The conversion table is found in the file vtag2arts.data.
%
% FORMAT:   name = vtag2arts( tag )
%
% OUT:      name    species name according to ARTS
% IN:       tag     Verdandi/HITRAN tag number
%------------------------------------------------------------------------

% HISTORY: 1996.11.26   Created by Patrick Eriksson (for Norns)
%          2001.03.26   Modified to AMI by Patrick Eriksson


function name = vtag2arts( tag )

% The function handled before several tags on the same time,
% and parts of this functionality is left in the function.
% Just if you wonder ;-)

niso	= 5;	%Maximum number of isotopes
smax	= 15;	%Maximum length of a name string


iso = rem( tag, 10 );
mo  = round( (tag-iso)/10 ); 


if iso > niso
  error('You selected a too high isotope number.');
end


fid	= fopen('vtag2arts.data');
if fid<0
  fprintf('tags2species:\n');
  fprintf('Could not open the file species.name');
  return
end

S	= 32*ones(length(mo),smax);
S	= setstr(S);

s	= 'a';
i	= 0;
n	= (mo-1)*(niso+1)+iso;
nmax	= max(n);
maxl	= -1;

while (s(1)~=-1) & (i<nmax)
  s = fgets(fid);
  i = i+1;
  ni = find(i==n);
  if ~isempty(ni)
    for j = 1:length(ni)
      ls = length(s);
      S(ni(j),1:ls) = s;
      if maxl<ls, maxl=ls; end
    end
  end
end

fclose(fid);

if maxl < 1
  error('The selected tag is not a valid option.');
end

name	= S(:,1:(maxl-1));
