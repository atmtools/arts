%-----------------------------------------------------------------------------
% NAME:     answer_is_yes
%
%           Prompts a question and only allows 'y' or 'n' as answer.
%
%           The question string S shall only contain the question without any
%           ? character. If s='Do you like money', the prompted question is:
%
%              Do you like money (y/n)?: 
%
% FORMAT:   bool = answer_is_yes( s )
%
% OUT:      bool   If anser is y, bool=1, else bool=0.
% IN:       s      String with question.
%-----------------------------------------------------------------------------

% HISTORY: 2002-03-10  Created by Patrick Eriksson


function bool = answer_is_yes( s )

response = 'w';

while length(response)~=1 | (response~='y' & response~='n')
  response = lower( input([s,' (y/n)?: '],'s') );
end

if response == 'y'
  bool = 1;
else
 bool = 0;
end


