function [reciprocalXgrid, reciprocalYgrid, reciprocalZgrid] = reciprocal_lab_space_coordinates(S, shape)
fprintf('reciprocal_lab_space_coordinates\n...creating orthogonal reciprocal lab space coordinates...');
% orthogonal meshgrid coordinates for reciprocal lab space
[x_pix, y_pix, z_pix] = size(shape);
[N1grid, N2grid, N3grid] = meshgrid(-(x_pix-1)/2:(x_pix-1)/2, -(y_pix-1)/2:(y_pix-1)/2, -(z_pix-1)/2:(z_pix-1)/2);

reciprocalXgrid = 2*pi()*N1grid/S.lambda*S.d/S.D; % equal to 2*pi/S.p_sam/S.N                     
reciprocalYgrid = 2*pi()*N2grid/S.lambda*S.d/S.D;
reciprocalZgrid = 2*pi()*N3grid/S.lambda*S.d/S.D;
fprintf('\n...done\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end