function [q_x, q_y, q_z] = RLS_basis_vectors(S)
% making RLS basis vectors
fprintf('RLS_basis_vectors\n...making RLS basis vectors...');
% making q_x, q_y, q_z lab space vectors
q_x = 2*pi/S.lambda*S.d/S.D*[1; 0; 0];
q_y = 2*pi/S.lambda*S.d/S.D*[0; 1; 0];
q_z = 2*pi/S.lambda*S.d/S.D*[0; 0; 1];
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end