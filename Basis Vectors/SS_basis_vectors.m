function [x_sam, y_sam, z_sam] = SS_basis_vectors(S)
% making SS basis vectors
fprintf('SS_basis_vectors\n...making SS basis vectors...');
% making x_sam , y_sam , z_sam sample space vectors
x_sam = S.R_xyz*S.p_sam*[1; 0; 0];
y_sam = S.R_xyz*S.p_sam*[0; 1; 0];
z_sam = S.R_xyz*S.p_sam*[0; 0; 1];
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end