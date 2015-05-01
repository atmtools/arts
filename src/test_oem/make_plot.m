% Script to plot benchmark results, assumes filename and
% plot_title to be in the workspace.

[n, t, t_m] = textread( filename, '', 'commentstyle', 'shell' );

f = figure( 'visible', 'off' )
set( gca, 'visible', 'off' )

plot( n, t, 'b-o', n, t_m, 'r-o' );
xlabel( 'Matrix Size' );
ylabel( 'CPU Time [ms]');
title( plot_title );

[path, name, ext] = fileparts( filename );
saveas( gcf, [name,'.png'] );