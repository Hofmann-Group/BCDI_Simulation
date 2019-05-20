function [N1gridp, N2gridp, N3gridp] = non_orthogonal_coordinates(S)
% makes non-orthogonal meshgrid coordinates based on the orthogonal coordinates and simulation geometry
fprintf('\nnon_orthogonal_coordinates');
%% Map from lab recoprocal space to detector reciprocal space:
fprintf('\n...creating non-orthogonal coordinates...');
% making q_x, q_y, q_z reciprocal lab space basis vectors
q_x = 2*pi()/S.lambda*S.d/S.D*[1; 0; 0];
q_y = 2*pi()/S.lambda*S.d/S.D*[0; 1; 0];
q_z = 2*pi()/S.lambda*S.d/S.D*[0; 0; 1];

% making q'_1, q'_2, q'_3 detector reciprocal space vectors
q_1p = S.R_dqp_12*2*pi()/S.lambda*S.d/S.D*[1; 0; 0];
q_2p = S.R_dqp_12*2*pi()/S.lambda*S.d/S.D*[0; 1; 0];
q_3p = S.R_dqp_3*S.Q_lab-S.Q_lab;

% make the T_RLS_to_DRS matrix
T_RLS_to_DRS = [q_x, q_y, q_z]\[q_1p, q_2p, q_3p]; % equivalent to [a_1, a_2, a_3]*inv(b_1, b_2, b_3])

% map orthogonal coordinates to non-orthogonal coordinates
N1gridp = T_RLS_to_DRS(1,1)*S.N1grid + T_RLS_to_DRS(1,2)*S.N2grid + T_RLS_to_DRS(1,3)*S.N3grid;
N2gridp = T_RLS_to_DRS(2,1)*S.N1grid + T_RLS_to_DRS(2,2)*S.N2grid + T_RLS_to_DRS(2,3)*S.N3grid;
N3gridp = T_RLS_to_DRS(3,1)*S.N1grid + T_RLS_to_DRS(3,2)*S.N2grid + T_RLS_to_DRS(3,3)*S.N3grid;
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end