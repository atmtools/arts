function fig_absbasis(ztan,zp1,zp2,zp3)
% Give altitudes in km


z   = -1:0.01:1;

n   = (length(z)-1)/2;
dy  = 1/n;

y   = [0:dy:1,(1-dy):-dy:0];


figure(1),clf
plot(y,zp1+z),hold on
plot(y,zp2+z,'--')
plot(y,zp3+z,'-.')
xlabel2('Basis function',14,1);
ylabel2('Altitude [km]',14,1);
axis([0 1.1 ztan-1 zp3+1.5])


figure(2),clf
plot(z2l(ztan*1e3,(zp1+z)*1e3)/1e3,y),hold on
plot(z2l(ztan*1e3,(zp2+z)*1e3)/1e3,y,'--')
plot(z2l(ztan*1e3,(zp3+z)*1e3)/1e3,y,'-.')
ylabel2('Basis function',14,1);
xlabel2('Distance from z_{tan} [km]',14,1);
axis([0 z2l(ztan*1e3,(zp3+max(z))*1e3)/1e3+10 0 1.1])


%figure(1),print fig_absbasis_z.eps -deps
%figure(2),print fig_absbasis_l.eps -deps



function l = z2l(ztan,z)

R = 6837e3;

l = sqrt((R+z).^2-(R+ztan)^2);
