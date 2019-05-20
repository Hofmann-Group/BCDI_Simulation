function [xp, yp, zp] = DCS_basis_vectors(S)
% making DCS basis vectors
fprintf('DCS_basis_vectors\n...making DCS basis vectors...');
try
    [N1, N2, N3] = size(S.SS_shape);
catch
    fprintf('\n...size of array set to N1 = N2 = N3 = 256...');
    N1 = 256;
    N2 = 256;
    N3 = 256;
end

% making q_1', q_2', q_3' detector reciprocal space vectors
q_1p = 2*pi/S.lambda*S.d/S.D*S.R_dqp_12*[1; 0; 0];
q_2p = 2*pi/S.lambda*S.d/S.D*S.R_dqp_12*[0; 1; 0];
q_3p = S.R_dqp_3*S.Q_lab-S.Q_lab;

% make x', y', and z' detector conjugated space vectors, adapted from Berenguer et al. PRB 88, 144101 (2013).
V_DRS = dot(cross(q_3p, q_2p), q_1p)*N1*N2*N3;
xp = 2*pi*cross(N2*q_2p, N3*q_3p)./V_DRS;
yp = 2*pi*cross(N3*q_3p, N1*q_1p)./V_DRS;
zp = 2*pi*cross(N1*q_1p, N2*q_2p)./V_DRS;
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end