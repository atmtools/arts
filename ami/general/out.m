%------------------------------------------------------------------------
% NAME:    out
%
%          Prints messages. Mimics the ARTS output streams.
%          The verbosity is determined by the global variable 
%          REPORT_LEVEL. 
%
%          The function puts a frame around the text. To begin a frame
%          put s=1, to end s=-1. A braking line is obtained by s=0.
%
%          An example. The command
%
%            out(1,1);out(1,'Heading');out(1,0);out(2,'Some text');out(1,-1)
%
%          gives
%
%           /-----------------------------------------------------------------\
%           | Heading                                                         |
%           |-----------------------------------------------------------------|
%           |   Some text                                                     |
%           \-----------------------------------------------------------------/
%
%          The function can also be used to determine if figures etc.
%          shall be produced (the last format example).
%
% FORMAT:  out( level, s )
%            or
%          do_output = out( level )
%
% OUT:     do_output   Boolean. True if level <= REPORT_LEVEL
% IN:      level       Report level for the message.
%          s           Message. Can be a string matrix (see str2mat).
%                      If s is an integer, different vertical lines are
%                      produced (see above).
%------------------------------------------------------------------------

% HISTORY: 00.12.24  Created by Patrick Eriksson. 


function do_output = out( level, s )


global REPORT_LEVEL


ncols = 70;


if ~isempty(REPORT_LEVEL)  & (level<=REPORT_LEVEL)
  do_output = 1;
else
  do_output = 0;
  return
end


if nargin == 1
  return;
end

if ischar(s)

  for i = 1:size(s,1)
  
    fprintf('| ');
  
    %Indention
    for j = 1:(level-1)
      fprintf('  ');
    end
  
    fprintf('%s',s(i,:));
  
    for j = 1:(ncols-length(s)-(level-1)*2-3)
      fprintf(' ');
    end
  
    fprintf('|\n');
  
  end


else

  %=== Start line
  if s > 0
    fprintf('\n/');
    for j = 1:(ncols-2)
      fprintf('-');
    end
    fprintf('\\\n');

  %=== Brake line
  elseif s == 0
    fprintf('|');
    for j = 1:(ncols-2)
      fprintf('-');
    end
    fprintf('|\n')

  %=== End line
  else
    fprintf('\\');
    for j = 1:(ncols-2)
      fprintf('-');
    end
    fprintf('/\n\n');
  end

end
