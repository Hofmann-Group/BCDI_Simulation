function [q_1p, q_2p, q_3p] = DRS_basis_vectors(S)
% making DRS basis vectors
fprintf('DRS_basis_vectors\n...making DRS basis vectors...');
% making q_1', q_2', q_3' detector reciprocal space vectors
q_1p = 2*pi/S.lambda*S.d/S.D*S.R_dqp_12*[1; 0; 0];
q_2p = 2*pi/S.lambda*S.d/S.D*S.R_dqp_12*[0; 1; 0];
q_3p = S.R_dqp_3*S.Q_lab-S.Q_lab;
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end