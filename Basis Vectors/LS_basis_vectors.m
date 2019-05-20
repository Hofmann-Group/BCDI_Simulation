function [x, y, z] = LS_basis_vectors(S)
% making LS basis vectors
fprintf('LS_basis_vectors\n...making LS basis vectors...');
% making x, y, z lab space vectors
x = S.p_sam*[1; 0; 0];
y = S.p_sam*[0; 1; 0];
z = S.p_sam*[0; 0; 1];
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end