function fig_corrfuns


x = -4:0.02:4;


figure(1),clf
plot(x,exp(-x.^2)),hold on
plot(x,exp(-abs(x)))
a = 1/(1-exp(-1));
plot([-a 0 a],[0 1 0])

ylabel2('Correlation',14,1);
xlabel2('Scaled frequency',14,1);

return

print fig_corrfuns.eps -deps
