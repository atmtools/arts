%------------------------------------------------------------------------
% NAME:    axes_nx1(n)
%
%          Sets tight axes for figures with a vertical row of subplots.
%
%          This function is useful when the horisontal grid of the subplots
%          is identical.
%
% RETURN:  -
% IN:      n       the number of sub-plots
%------------------------------------------------------------------------

% HISTORY: 2001-08-23  Copied from Urd (plotnx1) by Patrick Eriksson

function axes_nx1(n)

if n==1, return,end

x0	= 0.13;
xw	= 0.80;

if n==2
  yh	= 0.42;		%=== for n=2
  y1	= 0.10;
  y2	= 0.54;

elseif n==3
  yh	= 0.29;		%=== for n=3
  y1	= 0.10;
  y2	= 0.40;
  y3	= 0.70;

elseif n==4
  yh	= 0.20;		%=== for n=4
  y1	= 0.10;
  y2	= 0.32;
  y3	= 0.54;
  y4    = 0.76;

else
  error(['Values not defined for: ',int2str(n)]);
end


if n == 2
  subplot(n,1,2)
  set(gca,'Pos',[x0 y1 xw yh]);
  subplot(n,1,1)
  set(gca,'Pos',[x0 y2 xw yh]);
  set(gca,'XTickLabel','');

elseif n == 3
  subplot(n,1,3)
  set(gca,'Pos',[x0 y1 xw yh]);
  subplot(n,1,2)
  set(gca,'Pos',[x0 y2 xw yh]);
  set(gca,'XTickLabel','');
  subplot(n,1,1)
  set(gca,'Pos',[x0 y3 xw yh]);
  set(gca,'XTickLabel','');

elseif n == 4
  subplot(n,1,4)
  set(gca,'Pos',[x0 y1 xw yh]);
  subplot(n,1,3)
  set(gca,'Pos',[x0 y2 xw yh]);
  set(gca,'XTickLabel','');
  subplot(n,1,2)
  set(gca,'Pos',[x0 y3 xw yh]);
  set(gca,'XTickLabel','');
  subplot(n,1,1)
  set(gca,'Pos',[x0 y4 xw yh]);
  set(gca,'XTickLabel','');

end
