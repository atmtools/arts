%------------------------------------------------------------------------
% NAME:    plineleg
%
%          Makes a fixed legend for lines plotted with PLINE.
%          See also PLINE
%
% RETURN:  -
% IN:      s       string matrix, each row correspond to one line type
%          xstart  left end point of lines (norm. unit)
%          ystart  upper limit for lines (norm. unit)
%          --------(Below optional parameters)
%          xlength length of lines in legend (norm. unit)
%          tsize   font size
%          ystep   the vertical spacing between the line (norm. unit)
% USING:   -
%------------------------------------------------------------------------

% HISTORY: 990926  Created by Patrick Eriksson.


function plineleg(s,xstart,ystart,xlength,tsize,ystep)


if ~exist('xlength')
  xlength = 0.10;       % norm. length of lines in legend
end
if ~exist('tsize')
  tsize=10;             % size of text
end


%==========================
ax     = axis;
dy     = ax(4)-ax(3);
y1     = ax(3)+dy*ystart; 
dx     = ax(2)-ax(1);
x1     = ax(1)+dx*xstart;
x2     = x1+dx*xlength;
x3     = x2+dx*0.02;
if ~exist('ystep')
  ystep = tsize*dy*5/1000;
else
  ystep = dy*ystep;
end

n  = size(s,1);

for i = 1:n

  yp = y1 - (i-1)*ystep;
                                     % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  switch i                           % Changes of the line types shall also
    case 1                           % be done in PLINE
      plot([x1 x2],[yp yp],'-');
    case 2
      plot([x1 x2],[yp yp],'--');
    case 3
      plot([x1 x2],[yp yp],'-.');
    case 4 
      plot([x1 x2],[yp yp],':','LineWidth',1);
    case 5 
      plot([x1 x2],[yp yp],'-','LineWidth',1);
    case 6 
      plot([x1 x2],[yp yp],'--','LineWidth',1);
    case 7 
      plot([x1 x2],[yp yp],'-.','LineWidth',1);
    case 8 
      plot([x1 x2],[yp yp],':','LineWidth',2);
    case 9 
      plot([x1 x2],[yp yp],'-','LineWidth',2);
    case 10 
      plot([x1 x2],[yp yp],'--','LineWidth',2);
    case 11 
      plot([x1 x2],[yp yp],'-.','LineWidth',2);
    case 12
      plot([x1 x2],[yp yp],':','LineWidth',3);
    case 13
      plot([x1 x2],[yp yp],'-','LineWidth',3);
    case 14
      plot([x1 x2],[yp yp],'--','LineWidth',3);
    case 15
      plot([x1 x2],[yp yp],'-.','LineWidth',3);
    case 16
      plot([x1 x2],[yp yp],':','LineWidth',4);
    case 17
      plot([x1 x2],[yp yp],'-','LineWidth',4);
    case 18
      plot([x1 x2],[yp yp],'--','LineWidth',4);
    case 19
      plot([x1 x2],[yp yp],'-.','LineWidth',4);
    case 20
      plot([x1 x2],[yp yp],':','LineWidth',4);
    case 21
      plot([x1 x2],[yp yp],'-','LineWidth',4);
    case 22
      plot([x1 x2],[yp yp],'--','LineWidth',4);
    case 23
      plot([x1 x2],[yp yp],'-.','LineWidth',4);
    otherwise
      error('Unknown line number');
  end

    h=text2(x3,yp,s(i,:),tsize,0);

end
