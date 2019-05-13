function [realXgrid, realYgrid, realZgrid] = lab_space_coordinates(S, shape)
fprintf('lab_space_coordinates\n...creating orthogonal lab space coordinates...');
% orthogonal meshgrid coordinates for lab space                                  
[x_pix, y_pix, z_pix] = size(shape);
[N1grid, N2grid, N3grid] = meshgrid(-(x_pix-1)/2:(x_pix-1)/2, -(y_pix-1)/2:(y_pix-1)/2, -(z_pix-1)/2:(z_pix-1)/2);

realXgrid = S.p_sam*N1grid;
realYgrid = S.p_sam*N2grid;
realZgrid = S.p_sam*N3grid;
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end