%------------------------------------------------------------------------
% NAME:    pline
%
%          Plots lines with different line style. The line style is
%          determined by a sequential number.
%
% RETURN:  h       handle to the plot line
% IN:      x       x coordinates
%          y       y coordinates
%          i       sequential number
% USING:   -
%------------------------------------------------------------------------

% HISTORY: 990925  Created by Patrick Eriksson.
    

function h = pline(x,y,i)

switch i
  case 1
    h = plot(x,y,'-');
  case 2
    h = plot(x,y,'--');
  case 3
    h = plot(x,y,'-.');
  case 4 
    h = plot(x,y,':','LineWidth',1);
  case 5 
    h = plot(x,y,'-','LineWidth',1);
  case 6 
    h = plot(x,y,'--','LineWidth',1);
  case 7 
    h = plot(x,y,'-.','LineWidth',1);
  case 8 
    h = plot(x,y,':','LineWidth',2);
  case 9 
    h = plot(x,y,'-','LineWidth',2);
  case 10 
    h = plot(x,y,'--','LineWidth',2);
  case 11
    h = plot(x,y,'-.','LineWidth',2);
  case 12
    h = plot(x,y,':','LineWidth',3);
  case 13
    h = plot(x,y,'-','LineWidth',3);
  case 14 
    h = plot(x,y,'--','LineWidth',3);
  case 15
    h = plot(x,y,'-.','LineWidth',3);
  case 16
    h = plot(x,y,':','LineWidth',4);
  case 17
    h = plot(x,y,'-','LineWidth',4);
  case 18 
    h = plot(x,y,'--','LineWidth',4);
  case 19
    h = plot(x,y,'-.','LineWidth',4);
  case 20
    h = plot(x,y,':','LineWidth',5);
  case 21
    h = plot(x,y,'-','LineWidth',5);
  case 22 
    h = plot(x,y,'--','LineWidth',5);
  case 23
    h = plot(x,y,'-.','LineWidth',5);
  otherwise
    error('Unknown line number');
end
