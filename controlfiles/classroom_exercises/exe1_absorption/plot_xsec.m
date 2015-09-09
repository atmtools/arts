%============================================================================
% NAME:    plot_xsec
%
%          this script reads and plots absorption cross sections 
%          corresponding to arts output from "absorption.arts"
%============================================================================

clear
close all


%% read in everything

abs_xsec_all = xmlLoad('./results/abs_xsec_per_species.xml');
abs_xsec = abs_xsec_all{1};
freq = xmlLoad('./results/f_grid.xml')/1e9;

species = xmlLoad('./results/species.xml');
mol_name = species{1};
temperature = xmlLoad('./results/temp.xml');
pressure = xmlLoad('./results/press.xml');


%% plot absorption cross section

plot(freq, abs_xsec, 'r', 'LineWidth',1);
xlabel('Frequency [ GHz ]', 'FontSize',14);
ylabel('Abs. cross section [ m^{2} ]', 'FontSize',14);
name = strrep(mol_name, '-*', '');
title(strcat( '\bf', name, ...
             ' P: ', num2str(pressure / 100.0), ' hPa ',...
             ' T: ', num2str(temperature), ' K'), 'FontSize',20);

%% save figure

pathname = './plots';
figname = strrep(strcat('plot_xsec_',mol_name,'_',num2str(pressure / 100.0),'hPa_',num2str(temperature),'K'), '-*','');
saveas(gcf,fullfile(pathname,figname{1}), 'epsc')