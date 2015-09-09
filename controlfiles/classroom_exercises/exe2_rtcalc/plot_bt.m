%==========================================================================
% NAME: plot_bt
% 
% 
%       this script plots zenith opacity and brightness temperature for different
%       species and different sensor positions. 
%       You need to run the arts control file "rtcalc.arts" first as this matlab
%       program will use the arts output variables. 
%==========================================================================

clear all
close all
clc

%% read input variables
bt = xmlLoad('./results/bt.xml'); 
freq = xmlLoad('./results/f_grid.xml');
mol_name  =  xmlLoad('./results/species.xml');
od = xmlLoad('./results/odepth_1D.xml'); 
sensor_pos = xmlLoad('./results/sensor_pos.xml'); 
sensor_los = xmlLoad('./results/sensor_los.xml'); 

%% figure of zenith opacity with logarithmic scale on y axis
% subplot(2,1,1)
semilogy(freq/10e8, od, 'Linewidth',1)
ylim([10e-6 100]);
xlabel('Frequency [ GHz ]', 'FontSize',12);
ylabel('zenith opacity', 'FontSize',12);
title('Zenith opacity','FontSize',14)
if length(mol_name)==3
text(10,3e-3,strrep(strcat('species: ',mol_name{1},', ',mol_name{2},' and ',mol_name{3}),'-*',''),'FontSize',10);
else
text(10,3e-3,strrep(strcat('species: ',mol_name{1},' and ',mol_name{2}),'-*',''),'FontSize',10);
end
text(10,5e-4,strcat('sensor position: ',num2str(sensor_pos/1000),'km'),'FontSize',10);
text(10,1e-4,strcat('sensor los: ',num2str(sensor_los),'degree'),'FontSize',10);

if length(mol_name) == 3
figname1 = strrep(strcat('opacity_',mol_name{1},'+',mol_name{2},'+',mol_name{3},'_',num2str(sensor_pos/1000),'km_',num2str(sensor_los),'deg'),'-*','');
else
figname1 = strrep(strcat('opacity_',mol_name{1},'+',mol_name{2},'_',num2str(sensor_pos/1000),'km_',num2str(sensor_los),'deg'),'-*','');  
end

path = './plots';
saveas(gcf, fullfile(path, figname1{1}), 'epsc')

%% figure of brithtness temperature for defined sensor position and line of
% sight (los)
% subplot(2,1,2)
figure()
plot(freq/10e8, bt, 'Linewidth',1, 'Color','r');
xlabel('Frequency [ GHz ]', 'FontSize',12);
ylabel('bt [ K ]', 'FontSize',12);
ylim([0 300]);
title('Brightness temperature','FontSize',14)
if length(mol_name)==3
text(10,100,strrep(strcat('species: ',mol_name{1},', ',mol_name{2},' and ',mol_name{3}),'-*',''),'FontSize',10);
else
text(10,100,strrep(strcat('species: ',mol_name{1},' and  ',mol_name{2}),'-*',''),'FontSize',10);     
end
text(10,130,strcat('sensor position: ',num2str(sensor_pos/1000),'km'),'FontSize',10);
text(10,160,strcat('sensor los: ',num2str(sensor_los),'degree'),'FontSize',10);

if length(mol_name) == 3
figname2 = strrep(strcat('bt_',mol_name{1},'+',mol_name{2},'+',mol_name{3},'_',num2str(sensor_pos/1000),'km_',num2str(sensor_los),'deg'),'-*','');
else
figname2 = strrep(strcat('bt_',mol_name{1},'+',mol_name{2},'_',num2str(sensor_pos/1000),'km_',num2str(sensor_los),'deg'),'-*','');  
end
%% save figure

saveas(gcf, fullfile(path, figname2{1}), 'epsc')


 