%------------------------------------------------------------------------
% NAME:   subfolders
%
%         Returns the name of the folders inside a folder.
%
%         The following folders are not returned:
%            .*  (any folder starting with a dot)
%            CVS
%
% FORMAT: names = subfolders( mainfolder )
%
% OUT:    names        Name on the subfolders as cell array of strings.
% IN:     mainfolder   Name in the folder where to look.
%------------------------------------------------------------------------

% HISTORY: 2002.01.08  Created by Patrick Eriksson.


function names = subfolders( mainfolder )


D = dir( mainfolder );

names = cell( length(D), 1 );

n = 0;


for i = 1:length(D)

  if D(i).isdir

    name = D(i).name;

    if ~strncmp(name,'.',1) & ~strcmp(name,'CVS')

      n = n + 1;

      names{n} = name;

    end

  end

end


names = names(1:n);