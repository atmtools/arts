function k = amenu(xHeader,varargin);
%MENU   Generate an ASCII menu of choices for user input.
%
%   CHOICE = MENU(HEADER, ITEM1, ITEM2, ... ) displays the HEADER
%   string followed in sequence by the menu-item strings: ITEM1, ITEM2,
%   ... ITEMn. Returns the number of the selected menu-item as CHOICE,
%   a scalar value. There is no limit to the number of menu items.
%
%   CHOICE = MENU(HEADER, ITEMLIST) where ITEMLIST is a string cell
%   array is also a valid syntax.
%
%   Command window example:
%   >> K = MENU('Choose a color','Red','Blue','Green')
%   displays on the screen:
%
%     ----- Choose a color -----
%
%        1) Red
%        2) Blue
%        3) Green
%
%        Select a menu number:
%
%   The number entered by the user in response to the prompt is
%   returned as K (i.e. K = 2 implies that the user selected Blue).

% This function is an excerpt from the Matlab function menu, where only
% the ASCII version is left.
%
% Patrick Eriksson 2001.08.07


%=========================================================================
% Check input
%-------------------------------------------------------------------------
if nargin < 2,
    disp('MENU: No menu items to chose from.')
    k=0;
    return;
elseif nargin==2 & iscell(varargin{1}),
  ArgsIn = varargin{1}; % a cell array was passed in
else,
  ArgsIn = varargin;    % use the varargin cell array
end

%-------------------------------------------------------------------------
% Create the menu
%-------------------------------------------------------------------------
%
k = local_ASCIImenu( xHeader, ArgsIn );

%#########################################################################
%   END   :  main function "menu"
%#########################################################################


%#########################################################################
%  BEGIN  :  local function local_ASCIImenu
%#########################################################################
function k = local_ASCIImenu( xHeader, xcItems )

% local function to display an ascii-generated menu and return the user's
% selection from that menu as an index into the xcItems cell array

%-------------------------------------------------------------------------
% Calculate the number of items in the menu
%-------------------------------------------------------------------------
numItems = length(xcItems);

%-------------------------------------------------------------------------
% Continuous loop to redisplay menu until the user makes a valid choice
%-------------------------------------------------------------------------
while 1,
    % Display the header
    disp(' ')
    disp(['----- ',xHeader,' -----'])
    disp(' ')
    % Display items in a numbered list
    for n = 1 : numItems
        disp( [ '      ' int2str(n) ') ' xcItems{n} ] )
    end
    disp(' ')
    % Prompt for user input
    k = input('Select a menu number: ');
    % Check input:
    % 1) make sure k has a value
    if isempty(k), k = -1; end;
    % 2) make sure the value of k is valid
    if  (k < 1) | (k > numItems) ...
        | ~strcmp(class(k),'double') ...
        | ~isreal(k) | (isnan(k)) | isinf(k),
        % Failed a key test. Ask question again
        disp(' ')
        disp('Selection out of range. Try again.')
    else
        % Passed all tests, exit loop and return k
        return
    end % if k...
end % while 1

%#########################################################################
%   END   :  local function local_ASCIImenu
%#########################################################################
