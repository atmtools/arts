function fig_calerror()

v    = [100e9,200e9,500e9];
t    = 1:300;
tcal = [2.7 300;80 300];

figure(1),clf
figure(2),clf
figure(3),clf
figure(4),clf

for k = 1:2
for j = 1:length(v)
  i    = planck(v(j),t);
  ical = planck(v(j),tcal(k,:));
  trj  = invrayjean(v(j),i);
  y    = tcal(k,1) + (tcal(k,2)-tcal(k,1))*(i-ical(1))/(ical(2)-ical(1));
  switch j
    case 1
      s = '-';
    case 2
      s = '--';
    case 3
      s = '-.';
  end
  figure((k-1)*2+1)
  plot(t,t-y,s),hold on
  figure((k-1)*2+2)
  plot(t,trj-y,s),hold on
end
end

for k=1:2
figure((k-1)*2+1)
xlabel('Calibrated temperature (T_{cal})','FontSize',12)
ylabel('T_b - T_{cal}  [K]','FontSize',12)
legend('100 GHz','200 GHz','500 GHz',4)
%axis([0 300 -10 0.5])
if k==1
  print fig_calerror_tb_1.eps -deps
else
  print fig_calerror_tb_2.eps -deps
end
end


for k=1:2
figure((k-1)*2+2)
xlabel('Calibrated temperature (T_{cal})','FontSize',12)
ylabel('T_{rj} - T_{cal}  [K]','FontSize',12)
legend('100 GHz','200 GHz','500 GHz',4)
%axis([0 300 -12 0])
if k==1
  print fig_calerror_trj_1.eps -deps
else
  print fig_calerror_trj_2.eps -deps
end
end


