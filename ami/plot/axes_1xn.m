%------------------------------------------------------------------------
% NAME:    axes_1xn(n)
%
%          Sets tight axes for figures with a horisontal row of subplots.
%
%          This function is useful when the vertical grid of the subplots
%          is identical.
%
% RETURN:  -
% IN:      n       the number of sub-plots
%------------------------------------------------------------------------

% HISTORY: 2001-08-23  Copied from Urd (plot1xn) by Patrick Eriksson

function axes_1xn(n)

if n==1, return,end

y0	= 0.1;
yh	= 0.82;

if n==2
  w	= 0.45;		%=== for n=2
  x1	= 0.08;
  x2	= 0.54;

elseif n==3
  w	= 0.30;		%=== for n=3
  x1	= 0.07;
  x2	= 0.38;
  x3	= 0.69;

elseif n==4
  w	= 0.22;		%=== for n=4
  x1	= 0.08;
  x2	= 0.31;
  x3	= 0.54;
  x4	= 0.77;

else
  error(['Values not defined for: ',int2str(n)]);
end


if n == 2
  subplot(1,n,1)
  set(gca,'Pos',[x1 y0 w yh]);
  subplot(1,n,2)
  set(gca,'Pos',[x2 y0 w yh]);
  set(gca,'YTickLabel','');

elseif n == 3
  subplot(1,n,1)
  set(gca,'Pos',[x1 y0 w yh]);
  subplot(1,n,2)
  set(gca,'Pos',[x2 y0 w yh]);
  set(gca,'YTickLabel','');
  subplot(1,n,3)
  set(gca,'Pos',[x3 y0 w yh]);
  set(gca,'YTickLabel','');

elseif n == 4
  subplot(1,n,4)
  set(gca,'Pos',[x4 y0 w yh]);
  set(gca,'YTickLabel','');
  subplot(1,n,3)
  set(gca,'Pos',[x3 y0 w yh]);
  set(gca,'YTickLabel','');
  subplot(1,n,2)
  set(gca,'Pos',[x2 y0 w yh]);
  set(gca,'YTickLabel','');
  subplot(1,n,1)
  set(gca,'Pos',[x1 y0 w yh]);

end
