%==========================================================================
% NAME: plot_olr
%
%==========================================================================

clear
close all
%% read input data

intensity = xmlLoad('./results/olr.xml');
freq = xmlLoad('./results/f_grid.xml');
t_surf = xmlLoad('./results/t_surface.xml');
scene = xmlLoad('./results/scene.xml');

wave_num   =  freq / 30e9;

%% calculate integrated spectra

int_spectra  =  trapz( freq, intensity );

%% Caculate Planck curves for different BB temperatures

bb_inten_Ts   =  planck( freq, t_surf);
bb_inten_275  =  planck( freq, 275.0);
bb_inten_250  =  planck( freq, 250.0);
bb_inten_225  =  planck( freq, 225.0);

%% Plotting...

set (gca,'Fontsize',15)

plot(wave_num, intensity);
xlabel('Wave number [ cm^{-1} ]');
ylabel('Intensity [ Wm^{-2}Hz^{-1}Sr^{-1} ]');

hold on

plot(wave_num, bb_inten_Ts,  'k');
plot(wave_num, bb_inten_275, 'r');
plot(wave_num, bb_inten_250, 'g');
plot(wave_num, bb_inten_225, 'm');

tg_str  =  num2str( t_surf, '%5.1f');

 legend( 'Intensity', strcat(tg_str, ' K (T_s)'), '275 K', '250 K', '225 K' )
legend( 'boxoff' )

text(1500, 3e-12, strcat( 'OLR = ',  num2str(int_spectra, '%5.2f'), ' Wm^{-2}Sr^{-1}'));

%% Save figure 

pathname = './plots';
figname = ['olr' scene '.eps'];
saveas (gcf, fullfile(pathname, figname), 'epsc');




