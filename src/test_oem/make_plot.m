% Script to plot benchmark results, assumes filename and
% plot_title to be in the workspace.

[n, t_blas, t_arts, t_m] = textread( filename, '', 'commentstyle', 'shell' );

f = figure( 'visible', 'off' )
set( gca, 'visible', 'off' )
fax = gca;

plot( fax, n, t_blas, 'b-o', n, t_arts , 'r-o', n, t_m, 'g-o' );
xlabel( 'Matrix Size' );
ylabel( 'CPU Time [ms]');
legend( 'BLAS', 'arts', 'Matlab')
title( plot_title );

[path, name, ext] = fileparts( filename );
saveas( gcf, [name,'.png'] );