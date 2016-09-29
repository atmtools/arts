%==========================================================================
%NAME: plot_jacobian
%
%      this script uses the output variables of the arts controlfile
%      "jacobian.arts". it plots nadir brightness temperature at the top of 
%      the atmosphere, zenith opacity, the water vapour jacobian and the 
%      opacity between altitude z and the top of the atmosphere. Different 
%      frequencies can be selected by changing the variable "freq_ind" into 
%      a number between 1 and 110.
%
%==========================================================================      

clear;
close all;

%% select frequency
freq_ind=10;

%% Read in everything
freq = xmlLoad('./results/f_grid.xml');
iy_aux = xmlLoad('./results/iy_aux.xml');
opacity = xmlLoad('./results/y_aux.xml');
bt = xmlLoad('./results/bt.xml');
jac = xmlLoad('./results/jacobian.xml');
alt = xmlLoad('./results/z_field.xml');
mol_name = xmlLoad('./results/species.xml');

press = iy_aux{1};
abs = iy_aux{2};
opacity = opacity{1};

%% Calculate the level where opacity becomes one

nlevls   =  length(alt);
nfreqs   =  length(freq);

%% Flipping the altitude dimension so that I can start the 
%% integration from TOA

abs     =  flipdim(abs, 2);
z_fine  =  flipdim(alt, 1);

z_chapmann  =  zeros(nfreqs, 1);

for ifreq = 1 : nfreqs
  tau  =  0.0;
  abs_per_freq  =  squeeze(abs(ifreq, :));
  
  for ilev = 1 : nlevls - 1
    abs_layer  =  0.5 * ( abs_per_freq(ilev) + abs_per_freq(ilev + 1));
    tau  =  tau + abs_layer * (z_fine(ilev) - (z_fine(ilev + 1)));
      
    if tau > 1
      z_chapmann(ifreq, 1)  =  z_fine(ilev) / 1000.0;
      break
    end
  
  end
end

%% Plotting ...
 
if (freq_ind == 0)
  disp('');
  disp('Error!!!: freq_ind must be either -1 or > 0 !!!');
  disp('Anyway, just see the plot of heigts for tau = 1');
  disp('');

  figure(1);
  plot(freq/10^9, z_chapmann)
  title('z(\tau = 1)')
  xlabel('Frequency [ GHz ]')
  ylabel('Altitude [ km ]')

  return
end
  
figure(1)

  subplot(2, 2, 1);
  plot(freq/10^9, bt)
  xlabel('\bf Frequency [ GHz ]');
  ylabel('\bf T_B [ K ]');
  title('\bf Brightness Temperature')

  if freq_ind > 0
    hold on 
    plot(freq(freq_ind)/10^9, bt(freq_ind), 'or')
    text( 152, 280, strcat('\bf', num2str(freq(freq_ind)/10^9, ...
                                          '%6.2f'), ' \bf', 'GHz'))
    hold off
  end
  
  subplot(2, 2, 2);
  semilogy(freq/10^9, opacity)
  xlabel('\bf Frequency [ GHz ]');
  ylabel('\bf Zenith Opacity [  ]');
  title('\bf Zenith Opacity')

  if freq_ind > 0
    hold on 
    plot(freq(freq_ind)/10^9, opacity(freq_ind), 'or')
    hold off
  end

if freq_ind > 0
  subplot(2, 2, 3);
  plot(jac(freq_ind, :), alt / 1000.0) 
  xlim([-1 7])
  ylim([0 70])
  xlabel('\bf Jacobian [ K / 1 ]');
  ylabel('\bf Altitude [ km ]');
  title('\bf Water Vapor Jacobian')

%% Calculate opacity between each level and TOA

  opac_layer_space  =  zeros(nlevls, 1);

  abs_per_freq  =  squeeze(abs(freq_ind, :));
  
  tau  = 0;
  for ilev = 1 : nlevls - 1
    abs_layer  =  0.5 * ( abs_per_freq(ilev) + abs_per_freq(ilev + 1) );

    tau  =  tau + abs_layer * (z_fine(ilev) - (z_fine(ilev + 1)));
    
    opac_layer_space(ilev + 1, 1)  =  tau;
  end

  subplot(2, 2, 4);
  semilogx(opac_layer_space, z_fine / 1000.0)
  xlim([10^-13 10^3])
  set(gca,'XTick',[10.^-13 10.^-10 10.^-5 10.^0 10.^2]);
  ylim([0 70]) 
  xlabel('\bf Opacity_{layer--space} [  ]');
  ylabel('\bf Altitude [ km ]');
  title('\bf Opacity from a layer to space')
  
  hold on
  plot([1, 1], [0, 20], '--')
end  

%% save figure

pathname = './plots';
if freq_ind < 0
    figname = 'bt_op_part1.eps';
    saveas (gcf, fullfile(pathname, figname), 'epsc');
else
    figname = strcat('jac_',num2str(round((freq(freq_ind)/10^9))),'GHz.eps');
    saveas (gcf, fullfile(pathname, figname), 'epsc');
end




