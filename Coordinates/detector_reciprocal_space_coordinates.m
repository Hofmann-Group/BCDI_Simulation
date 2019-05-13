function [det_reciprocalXgrid, det_reciprocalYgrid, det_reciprocalZgrid] = detector_reciprocal_space_coordinates(S)
fprintf('\ndetector_reciprocal_space_coordinates');
%% Declaring variables from shape structure
% detector parameters
D = S.D;
d = S.d;

% wavelength
lambda = S.lambda;

% reciprocal lattice grids
SrecipXgrid = S.reciprocalXgrid;
SrecipYgrid = S.reciprocalYgrid;
SrecipZgrid = S.reciprocalZgrid;

% rotation matrices
R_dqp_12 = S.R_dqp_12;
R_dqp_3 = S.R_dqp_3;

% scattering vector
Q_lab = S.Q_lab;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map from lab recoprocal space to detector reciprocal space:
fprintf('\n...creating detector reciprocal space coordinates...');
% making q_x, q_y, q_z reciprocal lab space basis vectors
q_x = 2*pi()/lambda*d/D*[1; 0; 0];
q_y = 2*pi()/lambda*d/D*[0; 1; 0];
q_z = 2*pi()/lambda*d/D*[0; 0; 1];

% making q'_1, q'_2, q'_3 detector reciprocal space vectors
q_1p = R_dqp_12*2*pi()/lambda*d/D*[1; 0; 0];
q_2p = R_dqp_12*2*pi()/lambda*d/D*[0; 1; 0];
q_3p = (R_dqp_3*Q_lab-Q_lab);

% make the T_RLS_to_DRS matrix
T_RLS_to_DRS = [q_x, q_y, q_z]\[q_1p, q_2p, q_3p]; % equivalent to [a_1, a_2, a_3]*inv(b_1, b_2, b_3])

% map orthogonal coordinates to non-orthogonal coordinates
det_reciprocalXgrid = T_RLS_to_DRS(1,1)*SrecipXgrid + T_RLS_to_DRS(1,2)*SrecipYgrid + T_RLS_to_DRS(1,3)*SrecipZgrid;
det_reciprocalYgrid = T_RLS_to_DRS(2,1)*SrecipXgrid + T_RLS_to_DRS(2,2)*SrecipYgrid + T_RLS_to_DRS(2,3)*SrecipZgrid;
det_reciprocalZgrid = T_RLS_to_DRS(3,1)*SrecipXgrid + T_RLS_to_DRS(3,2)*SrecipYgrid + T_RLS_to_DRS(3,3)*SrecipZgrid;
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end