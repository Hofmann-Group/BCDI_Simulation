function [realXgrid, realYgrid, realZgrid] = orthogonal_coordinates(shape)
fprintf('orthogonal_coordinates\n...creating orthogonal coordinates...');
% orthogonal meshgrid coordinates                                 
[x_pix, y_pix, z_pix] = size(shape);
[N1grid, N2grid, N3grid] = meshgrid(-(x_pix-1)/2:(x_pix-1)/2, -(y_pix-1)/2:(y_pix-1)/2, -(z_pix-1)/2:(z_pix-1)/2);

realXgrid = N1grid;
realYgrid = N2grid;
realZgrid = N3grid;
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end