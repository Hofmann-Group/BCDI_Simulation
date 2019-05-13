function DCS_shape_LS = LS_to_DCS(S)
fprintf('LS_to_DCS\n...transforming LS shape to DCS shape...');
%% Declaring variables from shape structure
% lab space shape
LS_shape_SS = S.LS_shape_SS;

% detector parameters
D = S.D;
d = S.d;

% sample pixel size
p_sam = S.p_sam;

% wavelength
lambda = S.lambda;

% real space grid
realXgrid = S.realXgrid;
realYgrid = S.realYgrid;
realZgrid = S.realZgrid;

% rotation matrices
R_dqp_12 = S.R_dqp_12;
R_dqp_3 = S.R_dqp_3;

% scattering vector
Q_lab = S.Q_lab;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Making basis vectors
% finding size of the array and creating detector conjugated space dimensions
[N1, N2, N3] = size(LS_shape_SS);

% making x, y, z lab space vectors
x = p_sam*[1; 0; 0];
y = p_sam*[0; 1; 0];
z = p_sam*[0; 0; 1];

% making q'_1, q'_2, q'_3 detector reciprocal space vectors
q_1p = 2*pi/lambda*d/D*R_dqp_12*[1; 0; 0];
q_2p = 2*pi/lambda*d/D*R_dqp_12*[0; 1; 0];
q_3p = R_dqp_3*Q_lab-Q_lab;

% make x', y', and z' detector conjugated space basis vectors, adapted from Berenguer et al. PRB 88, 144101 (2013).
V_DRS = dot(cross(q_3p, q_2p), q_1p)*N1*N2*N3;
xp = 2*pi*cross(N2*q_2p, N3*q_3p)./V_DRS;
yp = 2*pi*cross(N3*q_3p, N1*q_1p)./V_DRS;
zp = 2*pi*cross(N1*q_1p, N2*q_2p)./V_DRS;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map from lab space to detector conjugated space:
% make the T_LS_to_DCS matrix
T_LS_to_DCS = [x, y, z]\[xp, yp, zp]; % equivalent to [a_1, a_2, a_3]*inv(b_1, b_2, b_3])

% map orthogonal coordinates to non-orthogonal coordinates
N1gridp = T_LS_to_DCS(1,1)*realXgrid + T_LS_to_DCS(1,2)*realYgrid + T_LS_to_DCS(1,3)*realZgrid;
N2gridp = T_LS_to_DCS(2,1)*realXgrid + T_LS_to_DCS(2,2)*realYgrid + T_LS_to_DCS(2,3)*realZgrid;
N3gridp = T_LS_to_DCS(3,1)*realXgrid + T_LS_to_DCS(3,2)*realYgrid + T_LS_to_DCS(3,3)*realZgrid;

% interpolate data in lab frame to detector conjugated frame
DCS_shape_LS = interp3(realXgrid, realYgrid, realZgrid, LS_shape_SS, N1gridp, N2gridp, N3gridp, 'linear',  0); % make any values outside master data zero.
fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end