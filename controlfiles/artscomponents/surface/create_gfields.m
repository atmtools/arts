% CREATE_GFIELDS
%
%   Creates a number of gridded fields, for providing input to
%   TestSurface.arts.

% 2012-06-25 Patrick Eriksson


function create_gfields
  

%- A very basic surface_scalar_reflectivity field:

% The test file sets lat_true to +54 deg, which gives r=0.8 if r set to zero
% at South pole and one at the North pole. This used for scalar case.
    
G.name      = 'Surface scalar reflectivity';
G.gridnames = { 'Frequency', 'Incidence angle', 'Latitude', 'Longitude' };
G.grids     = { 10e9, [0 90], [-90 90], [-360 360] };
G.dataname  = 'Scalar reflectivity';
r           = zeros( [1 1 2 2] );
r(:,:,2,:)  = 1;
G.data      = repmat( r, [1 length(G.grids{2}) 1 1] );
%
xmlStore( 'scalar_r_field_2angles.xml', G, 'GriddedField4' );



%- The same, but with frequencies and zenith angles triggering more complex
%  interpolations
  
G.grids{1}  = [1e9 20e9];
G.grids{2}  = 0:10:90;
G.data      = repmat( r, [2 length(G.grids{2}) 1 1] );
%
xmlStore( 'scalar_r_field_10angles.xml', G, 'GriddedField4' );


%- A very basic surface_reflectivity field:
%  (diagonal matrix with 0.8 everywhere)
  
G.name      = 'Surface reflectivity';
G.gridnames = { 'Frequency', 'Incidence angle', 'Latitude', 'Longitude',...
                'Stokes element', 'Stokes element' };
G.grids     = { {'1'}, [0 90], [-90 90], [-360 360], 1:4, 1:4 };
G.dataname  = 'Reflectivity';
r              = zeros( [1 1 1 1 4 4] );
r(1,1,1,1,:,:) = 0.8*eye(4);
G.data         = repmat( r, [1 2 2 2 1 1] );
%
xmlStore( 'r_field_2angles.xml', G, 'GriddedField6' );



%- A very basic complex_n field:
%  (data taht gibe [3,5,0) at +54 deg.
  
G.name      = 'Surface complex refractive index';
G.gridnames = { 'Frequency', 'Real/Imaginary', 'Latitude', 'Longitude' };
G.grids     = { 10e9, [1 2], [-90 90], [-360 360] };
G.dataname  = 'Complex refractivity';
r           = zeros( [1 2 2 2] );
r(:,1,1,:)  = 1;
r(:,1,2,:)  = 3.5;
G.data      = r;
%
xmlStore( 'complex_n_field.xml', G, 'GriddedField4' );
