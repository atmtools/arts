%-----------------------------------------------------------------------------
% NAME:     mkfig_1gridbox_2d
%
%           This is just a temporary function that will be removed from AMI
%           at a later stage.
%
% FORMAT:   mkfig_1gridbox_2d( r1, r2, r3, r4, alpha1, alpha3, rground1, 
%                                  rground2, bbg, r_p, alpha_p, psi_p, l_max )
%-----------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


function mkfig_1gridbox_2d( r1, r2, r3, r4, alpha1, alpha3, rground1, rground2, bbg, r_p, alpha_p, psi_p, l_max )


% Length unit is set to km
lscale = 1e3;


% Expand position vector is length is 1
n = max([ length(r_p), length(alpha_p), length(psi_p) ]);
%
if length(r_p) == 1
  r_p = r_p*ones(n,1);
end
if length(alpha_p) == 1
  alpha_p = alpha_p*ones(n,1);
end
if length(psi_p) == 1
  psi_p = psi_p*ones(n,1);
end


% Plot grid box
plot_1gridbox_2d( r1, r2, r3, r4, alpha1, alpha3, rground1, rground2 );



% Loop paths
for i = 1:length(r_p)

  [r,alpha,l_step,psi,r_istep,alpha_istep,r_tan,alpha_tan,ground,i_ground] = do_1gridbox_2d_geom( r1, r2, r3, r4, alpha1, alpha3, rground1, rground2, bbg, r_p(i), alpha_p(i), psi_p(i), l_max );

  fprintf('\n Path %d\n',i);
  fprintf('r_istep  = %d\n',r_istep);
  fprintf('al_istep = %d\n',alpha_istep);
  fprintf('l_step   = %.3f-%.3f km\n',min(l_step)/lscale,max(l_step)/lscale);
  fprintf('psi_q    = %.6f\n',last(psi));
  fprintf('ground   = %d\n',ground);
  fprintf('i_ground = %d\n',i_ground);
  if ~isnan( r_tan )
    fprintf('r_tan    = %.4f km\n',r_tan/1e3);
    fprintf('a_tan    = %.3f\n',alpha_tan);
  end

  % Start point
  [x1,y1] = cyl2cart( r_p(i), alpha_p(i) );
  plot( x1/lscale, y1/lscale, 'r*' );  

  for j = 1:length(r)
    [x2,y2] = cyl2cart( r(j), alpha(j) );
    plot( x2/lscale, y2/lscale, 'ro' );  
    plot( [x1,x2]/lscale, [y1,y2]/lscale, 'r-' );
    x1 = x2;
    y1 = y2;
  end
 
  if( ~isnan(r_tan) )
    [x1,y1] = cyl2cart( r_tan, alpha_tan );
    plot( x1/lscale, y1/lscale, 'r+' ); 
  end

end


ax = axis;
axis( [ax(1) ax(2)+0.1 ax(3) ax(4)+0.1 ] );
if (ax(2)-ax(1))<=5*(ax(4)-ax(3))
  axis equal
end
